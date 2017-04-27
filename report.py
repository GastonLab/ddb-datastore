#!/usr/bin/env python

import utils
import numpy as np
import sys
import getpass
import argparse
from collections import defaultdict

from toil.job import Job
from ddb import configuration
from ddb_ngsflow import pipeline
from variantstore import Variant
from variantstore import SampleVariant
from coveragestore import SampleCoverage
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider


def process_sample(job, config, sample, samples, addresses, authenticator, thresholds, callers):
    job.fileStore.logToMaster("Retrieving data for sample {}\n".format(sample))
    job.fileStore.logToMaster("Retrieving coverage data from database\n")
    connection.setup(addresses, "coveragestore", auth_provider=authenticator)

    report_data = dict()
    filtered_variant_data = defaultdict(list)
    off_target_amplicon_counts = defaultdict(int)
    target_amplicon_coverage = defaultdict(lambda: defaultdict(float))

    iterated = 0
    passing_variants = 0
    filtered_low_freq = 0
    filtered_low_depth = 0
    filtered_off_target = 0

    for library in samples[sample]:
        report_panel_path = "/mnt/shared-data/ddb-configs/disease_panels/{}/{}" \
                            "".format(samples[sample][library]['panel'], samples[sample][library]['report'])

        job.fileStore.logToMaster("{}: processing amplicons from file {}".format(library, report_panel_path))
        target_amplicons = utils.get_target_amplicons(report_panel_path)

        for amplicon in target_amplicons:
            coverage_data = SampleCoverage.objects.timeout(None).filter(
                SampleCoverage.sample == samples[sample][library]['sample_name'],
                SampleCoverage.amplicon == amplicon,
                SampleCoverage.run_id == samples[sample][library]['run_id'],
                SampleCoverage.library_name == samples[sample][library]['library_name'],
                SampleCoverage.program_name == "sambamba"
            )
            ordered_amplicons = coverage_data.order_by('amplicon', 'run_id').limit(coverage_data.count() + 1000)
            for result in ordered_amplicons:
                target_amplicon_coverage[amplicon]['num_reads'] = result.num_reads
                target_amplicon_coverage[amplicon]['mean_coverage'] = result.mean_coverage

        job.fileStore.logToMaster("{}: retrieving variants".format(library))
        variants = SampleVariant.objects.timeout(None).filter(
            SampleVariant.reference_genome == config['genome_version'],
            SampleVariant.sample == samples[sample][library]['sample_name'],
            SampleVariant.run_id == samples[sample][library]['run_id'],
            SampleVariant.library_name == samples[sample][library]['library_name'],
            SampleVariant.max_maf_all <= thresholds['max_maf']
        ).allow_filtering()

        num_var = variants.count()
        ordered = variants.order_by('library_name', 'chr', 'pos', 'ref', 'alt').limit(variants.count() + 1000)
        job.fileStore.logToMaster("{}: retrieved {} variants from database\n".format(library, num_var))
        job.fileStore.logToMaster("{}: classifying and filtering variants\n".format(library))

        for variant in ordered:
            iterated += 1
            if variant.amplicon_data['amplicon'] is 'None':
                filtered_off_target += 1
                off_target_amplicon_counts[variant.amplicon_data['amplicon']] += 1
            else:
                amplicons = variant.amplicon_data['amplicon'].split(',')
                assignable = 0
                for amplicon in amplicons:
                    if amplicon in target_amplicons:
                        assignable += 1
                        break
                if assignable:
                    match_variants = Variant.objects.timeout(None).filter(
                        Variant.reference_genome == config['genome_version'],
                        Variant.chr == variant.chr,
                        Variant.pos == variant.pos,
                        Variant.ref == variant.ref,
                        Variant.alt == variant.alt
                    ).allow_filtering()

                    num_matches = match_variants.count()
                    ordered_var = match_variants.order_by('sample', 'library_name', 'run_id').limit(num_matches + 1000)
                    vafs = list()
                    for var in ordered_var:
                        vaf = var.max_som_aaf
                        vafs.append(vaf)
                    variant.vaf_median = np.median(vafs)
                    variant.vaf_std_dev = np.std(vafs)

                    # Putting in to Tier1 based on COSMIC
                    if variant.cosmic_ids:
                        if variant.max_som_aaf < thresholds['min_saf']:
                            filtered_variant_data['tier1_fail_variants'].append(variant)
                            filtered_low_freq += 1
                        elif variant.max_depth < thresholds['depth']:
                            filtered_variant_data['tier1_fail_variants'].append(variant)
                            filtered_low_depth += 1
                        else:
                            filtered_variant_data['tier1_pass_variants'].append(variant)
                            passing_variants += 1
                        continue

                    # Putting in to Tier1 based on ClinVar not being None or Benign
                    if variant.clinvar_data['pathogenic'] != 'None':
                        if variant.clinvar_data['pathogenic'] != 'benign':
                            if variant.clinvar_data['pathogenic'] != 'likely-benign':
                                if variant.max_som_aaf < thresholds['min_saf']:
                                    filtered_variant_data['tier1_fail_variants'].append(variant)
                                    filtered_low_freq += 1
                                elif variant.max_depth < thresholds['depth']:
                                    filtered_variant_data['tier1_fail_variants'].append(variant)
                                    filtered_low_depth += 1
                                else:
                                    filtered_variant_data['tier1_pass_variants'].append(variant)
                                    passing_variants += 1
                                continue

                    if variant.severity == 'MED' or variant.severity == 'HIGH':
                        if variant.max_som_aaf < thresholds['min_saf']:
                            filtered_variant_data['tier3_fail_variants'].append(variant)
                            filtered_low_freq += 1
                        elif variant.max_depth < thresholds['depth']:
                            filtered_variant_data['tier3_fail_variants'].append(variant)
                            filtered_low_depth += 1
                        else:
                            filtered_variant_data['tier3_pass_variants'].append(variant)
                            passing_variants += 1
                        continue
                    else:
                        if variant.max_som_aaf < thresholds['min_saf']:
                            filtered_variant_data['tier4_fail_variants'].append(variant)
                            filtered_low_freq += 1
                        elif variant.max_depth < thresholds['depth']:
                            filtered_variant_data['tier4_fail_variants'].append(variant)
                            filtered_low_depth += 1
                        else:
                            filtered_variant_data['tier4_pass_variants'].append(variant)
                            passing_variants += 1
                        continue
                else:
                    filtered_off_target += 1
                    off_target_amplicon_counts[variant.amplicon_data['amplicon']] += 1

        job.fileStore.logToMaster("{}: iterated through {} variants\n".format(library, iterated))
        job.fileStore.logToMaster("{}: filtered {} off-target variants\n".format(library, filtered_off_target))
        job.fileStore.logToMaster("{}: filtered {} low-freq variants\n".format(library, filtered_low_freq))
        job.fileStore.logToMaster("{}: filtered {} low-depth variants\n".format(library, filtered_low_depth))

        job.fileStore.logToMaster("{}: passing {} tier 1 and 2 variants"
                                  "\n".format(library, len(filtered_variant_data['tier1_pass_variants'])))
        job.fileStore.logToMaster("{}: passing {} tier3 variants"
                                  "\n".format(library, len(filtered_variant_data['vus_pass_variants'])))
        job.fileStore.logToMaster("{}: passing {} tier 4 variants"
                                  "\n".format(library, len(filtered_variant_data['tier4_pass_variants'])))

    report_data['variants'] = filtered_variant_data
    report_data['coverage'] = target_amplicon_coverage

    return report_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)", default='report')
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    parser.add_argument('-d', '--min_depth', help='Minimum depth threshold for variant reporting', default=250.0)
    parser.add_argument('-t', '--min_somatic_var_freq', help='Minimum reportable somatic variant frequency',
                        default=0.01)
    parser.add_argument('-p', '--max_pop_freq', help='Maximum allowed population allele frequency', default=0.005)
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    args.logLevel = "INFO"

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    libraries = configuration.configure_samples(args.samples_file, config)

    samples = configuration.merge_library_configs_samples(libraries)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
    else:
        auth_provider = None

    thresholds = {'min_saf': args.min_somatic_var_freq,
                  'max_maf': args.max_pop_freq,
                  'depth': args.min_depth}

    callers = ("mutect", "platypus", "vardict", "scalpel", "freebayes", "pindel")

    sys.stdout.write("Processing samples\n")
    root_job = Job.wrapJobFn(pipeline.spawn_batch_jobs, cores=1)

    for sample in samples:
        sample_job = Job.wrapJobFn(process_sample, config, sample, samples, [args.address], auth_provider,
                                   thresholds, callers, cores=1)
        root_job.addChild(sample_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)

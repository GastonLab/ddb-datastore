#!/usr/bin/env python

import argparse
import getpass
import xlsxwriter
import sys
import csv

from cassandra.auth import PlainTextAuthProvider
from cassandra.cqlengine import connection
from ddb import configuration

import utils
from variantstore import SampleVariant
from coveragestore import SampleCoverage

from collections import defaultdict


def get_target_amplicons(filename):
    amplicons_list = list()
    with open(filename, "r") as bedfile:
        reader = csv.reader(bedfile, dialect='excel-tab')
        for row in reader:
            amplicons_list.append(row[3])

    return amplicons_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)")
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)

    args = parser.parse_args()

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    libraries = configuration.configure_samples(args.samples_file, config)

    samples = configuration.merge_library_configs_samples(libraries)
    preferred_transcripts = utils.get_preferred_transcripts(config['transcripts'])

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "variantstore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "variantstore")

    thresholds = {'min_saf': 0.01,
                  'max_maf': 0.005,
                  'depth': 200}

    callers = ("mutect", "platypus", "vardict", "scalpel", "freebayes", "pindel")

    sys.stdout.write("Processing samples\n")
    for sample in samples:
        sys.stdout.write("Processing variants for sample {}\n".format(sample))

        passing_variants = list()
        reportable_amplicons = list()

        filtered_no_amplicon = list()
        filtered_non_target_amplicon = list()
        filtered_no_requested_caller = list()

        target_amplicon_coverage = defaultdict(lambda: defaultdict(float))

        for library in samples[sample]:
            report_panel_path = "/mnt/shared-data/ddb-configs/disease_panels/{}/{}" \
                                "".format(samples[sample][library]['panel'], samples[sample][library]['report'])
            target_amplicons = get_target_amplicons(report_panel_path)

            sys.stdout.write("Processing variants for library {}\n".format(library))
            sys.stdout.write("Processing amplicons for library from file {}\n".format(report_panel_path))

            variants = SampleVariant.objects.timeout(None).filter(
                SampleVariant.reference_genome == config['genome_version'],
                SampleVariant.sample == samples[sample][library]['sample_name'],
                SampleVariant.run_id == samples[sample][library]['run_id'],
                SampleVariant.library_name == samples[sample][library]['library_name'],
                SampleVariant.max_som_aaf >= thresholds['min_saf'],
                SampleVariant.max_maf_all <= thresholds['max_maf'],
                ).allow_filtering()

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
                    reportable_amplicons.append(result)
                    target_amplicon_coverage[amplicon]['num_reads'] = result.num_reads
                    target_amplicon_coverage[amplicon]['mean_coverage'] = result.mean_coverage

            ordered_variants = variants.order_by('library_name', 'chr', 'pos',
                                                 'ref', 'alt').limit(variants.count() + 1000)

            for variant in ordered_variants:
                if variant.amplicon_data['amplicon']:
                    amplicons = variant.amplicon_data['amplicon'].split(',')
                    for amplicon in amplicons:
                        if amplicon in target_amplicons:
                            for caller in callers:
                                if caller in variant.callers:
                                    passing_variants.append(variant)
                                    break
                        else:
                            filtered_non_target_amplicon.append(variant)
                else:
                    filtered_no_amplicon.append(variant)

            sys.stdout.write("Retrieved {} total variants\n".format(variants.count()))

        sys.stdout.write("Sending {} variants to reporting (filtered {} variants for no amplicon data and {} for being"
                         " in a non-targeted amplicon)\n".format(len(passing_variants), len(filtered_no_amplicon),
                                                                 len(filtered_non_target_amplicon)))
        utils.write_sample_variant_report(args.report, sample, passing_variants, target_amplicon_coverage, callers,
                                          preferred_transcripts)

        sys.stdout.write("Writing coverage report\n")
        with open("{}_coverage_{}.txt".format(sample, args.report), "w") as coverage_report:
            coverage_report.write("Sample\tLibrary\tAmplicon\tNum Reads\tCoverage\n")
            for amplicon in reportable_amplicons:
                coverage_report.write("{}\t{}\t{}\t{}\t{}\n".format(amplicon.sample,
                                                                    amplicon.library_name,
                                                                    amplicon.amplicon,
                                                                    amplicon.num_reads,
                                                                    amplicon.mean_coverage))

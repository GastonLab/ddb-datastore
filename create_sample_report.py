#!/usr/bin/env python

import argparse
import getpass
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
    samples = configuration.configure_samples(args.samples_file, config)

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
        passing_variants = list()
        reportable_amplicons = list()
        filtered_no_amplicon = list()
        filtered_non_target_amplicon = list()
        filtered_no_requested_caller = list()

        target_amplicon_coverage = defaultdict(lambda: defaultdict(float))

        report_panel_path = "/mnt/shared-data/ddb-configs/disease_panels/{}/{}" \
                            "".format(samples[sample]['panel'], samples[sample]['report'])
        target_amplicons = get_target_amplicons(report_panel_path)

        sys.stdout.write("Running Cassandra query for sample {}\n".format(sample))

        variants = SampleVariant.objects.timeout(None).filter(
            SampleVariant.reference_genome == config['genome_version'],
            SampleVariant.sample == samples[sample]['sample_name'],
            SampleVariant.run_id == samples[sample]['run_id'],
            SampleVariant.library_name == samples[sample]['library_name'],
            SampleVariant.max_som_aaf >= thresholds['min_saf'],
            SampleVariant.max_maf_all <= thresholds['max_maf'],
            ).allow_filtering()

        for amplicon in target_amplicons:
            coverage_data = SampleCoverage.objects.timeout(None).filter(
                SampleCoverage.sample == samples[sample]['sample_name'],
                SampleCoverage.amplicon == amplicon,
                SampleCoverage.run_id == samples[sample]['run_id'],
                SampleCoverage.library_name == samples[sample]['library_name'],
                SampleCoverage.program_name == "sambamba"
            )

            for result in coverage_data:
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
        sys.stdout.write("Sending {} variants to reporting\n".format(len(passing_variants)))
        utils.write_sample_variant_report(args.report, sample, passing_variants, target_amplicon_coverage, callers)

        sys.stdout.write("Writing coverage report\n")
        with open("{}_coverage_{}.txt".format(sample, args.report), "w") as coverage_report:
            coverage_report.write("Sample\tLibrary\tAmplicon\tNum Reads\tCoverage\n")
            for amplicon in reportable_amplicons:
                coverage_report.write("{}\t{}\t{}\t{}\t{}\n".format(amplicon.sample,
                                                                    amplicon.library_name,
                                                                    amplicon.amplicon,
                                                                    amplicon.num_reads,
                                                                    amplicon.mean_coverage))

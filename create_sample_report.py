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
        target_amplicons = get_target_amplicons("/mnt/shared-data/ddb-configs/disease_panels/{}/{}"
                                                "".format(samples[sample]['panel'], samples[sample]['report']))
        target_amplicon_coverage = defaultdict(lambda: defaultdict(float))

        sys.stdout.write("Running Cassandra query for sample {}\n".format(sample))

        variants = SampleVariant.objects.timeout(None).filter(
            SampleVariant.reference_genome == config['genome_version'],
            SampleVariant.sample == samples[sample]['sample_name'],
            SampleVariant.run_id == samples[sample]['run_id'],
            SampleVariant.library_name == samples[sample]['library_name'],
            SampleVariant.max_som_aaf >= thresholds['min_saf'],
            SampleVariant.max_maf_all <= thresholds['max_maf'],
            SampleVariant.max_depth >= thresholds['depth']
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
                target_amplicon_coverage[amplicon]['num_reads'] = result.num_reads
                target_amplicon_coverage[amplicon]['mean_coverage'] = result.mean_coverage

        ordered_variants = variants.order_by('library_name', 'chr', 'pos',
                                             'ref', 'alt').limit(variants.count() + 1000)
        filtered_variants = list()
        for variant in ordered_variants:
            if variant.amplicon_data['amplicon']:
                amplicons = variant.amplicon_data['amplicon'].split(',')
                for amplicon in amplicons:
                    if amplicon in target_amplicons:
                        for caller in callers:
                            if caller in variant.callers:
                                filtered_variants.append(variant)
                                break

        sys.stdout.write("Retrieved {} total variants\n".format(variants.count()))
        sys.stdout.write("Sending {} variants to reporting\n".format(len(filtered_variants)))
        utils.write_sample_variant_report(args.report, sample, filtered_variants, target_amplicon_coverage, callers)

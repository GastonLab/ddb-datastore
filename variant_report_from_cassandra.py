#!/usr/bin/env python

import sys
import argparse
import utils
from ddb import configuration
from variantstore import SampleVariant
from cassandra.cqlengine import connection

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)")
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-v', '--variant_callers', help="Comma-delimited list of variant callers used")

    args = parser.parse_args()

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    connection.setup([args.address], "variantstore")

    thresholds = {'min_saf': 0.01,
                  'max_maf': 0.01,
                  'regions': config['actionable_regions']}

    callers = ['mutect', 'freebayes', 'scalpel', 'vardict', 'platypus', 'pindel']

    sys.stdout.write("Processing samples\n")
    for sample in samples:
        sys.stdout.write("Running Cassandra query for sample {}\n".format(sample))
        variants = SampleVariant.objects.timeout(None).filter(SampleVariant.sample == samples[sample]['sample_name'],
                                                SampleVariant.max_som_aaf >= thresholds['min_saf']).allow_filtering()
        ordered_variants = variants.order_by('library_name', 'reference_genome', 'chr', 'pos').limit(variants.count() + 1)
        sys.stdout.write("Retrieved {} total variants\n".format(variants.count()))
        sys.stdout.write("Running filters on sample variants\n")
        passing_variants = list()
        passed = 0
        iterated = 0
        for variant in ordered_variants:
            iterated += 1
            flag, info = utils.variant_filter(variant, callers, thresholds)
            if variant.max_som_aaf >= thresholds['min_saf']:
                passed += 1
                passing_variants.append((variant, flag, info))

        sys.stdout.write("Writing {} passing variants (of {}) to sample report\n".format(passed, iterated))
        utils.write_sample_variant_report(args.report, sample, passing_variants, args.variant_callers, thresholds)

#!/usr/bin/env python

import argparse
import getpass
import sys

from cassandra.auth import PlainTextAuthProvider
from cassandra.cqlengine import connection
from ddb import configuration

from scripts import utils
from variantstore import SampleVariant

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)")
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    parser.add_argument('-v', '--variant_callers', help="Comma-delimited list of variant callers used")

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

    thresholds = {'min_saf': 0.02,
                  'max_maf': 0.005,
                  'depth': 500,
                  'regions': config['actionable_regions']}

    callers = ['mutect', 'freebayes', 'scalpel', 'vardict', 'platypus', 'pindel']

    sys.stdout.write("Processing samples\n")
    for sample in samples:
        sys.stdout.write("Running Cassandra query for sample {}\n".format(sample))
        variants1 = SampleVariant.objects.timeout(None).filter(SampleVariant.sample == samples[sample]['sample_name'],
                                                               SampleVariant.run_id == samples[sample]['run_id'],
                                                               SampleVariant.library_name == samples[sample]['library1'],
                                                               SampleVariant.max_som_aaf >= thresholds['min_saf'],
                                                               SampleVariant.max_maf_all <= thresholds['max_maf'],
                                                               SampleVariant.max_depth >= thresholds['depth']
                                                               ).allow_filtering()

        variants2 = SampleVariant.objects.timeout(None).filter(SampleVariant.sample == samples[sample]['sample_name'],
                                                               SampleVariant.run_id == samples[sample]['run_id'],
                                                               SampleVariant.library_name == samples[sample]['library1'],
                                                               SampleVariant.max_som_aaf >= thresholds['min_saf'],
                                                               SampleVariant.max_maf_all <= thresholds['max_maf'],
                                                               SampleVariant.max_depth >= thresholds['depth']
                                                               ).allow_filtering()

        ordered_variants1 = variants1.order_by('library_name', 'reference_genome',
                                               'chr', 'pos').limit(variants1.count() + 1000)

        ordered_variants2 = variants2.order_by('library_name', 'reference_genome',
                                               'chr', 'pos').limit(variants2.count() + 1000)

        sys.stdout.write("Retrieved {} total variants from library {}\n".format(variants1.count(),
                                                                                samples[sample]['library1']))
        sys.stdout.write("Retrieved {} total variants from library {}\n".format(variants2.count(),
                                                                                samples[sample]['library2']))

        sys.stdout.write("Running filters on sample variants\n")
        dual_lib_variants = list()
        passing_variants = list()
        iterated = 0

        for variant in ordered_variants1:
            iterated += 1
            flag, info = utils.variant_filter(variant, callers, thresholds)
            if info['dual']:
                dual_lib_variants.append((variant, flag, info))
                passing_variants.append((variant, flag, info))
            else:
                passing_variants.append((variant, flag, info))

        for variant in ordered_variants2:
            iterated += 1
            flag, info = utils.variant_filter(variant, callers, thresholds)
            if info['dual']:
                dual_lib_variants.append((variant, flag, info))
            else:
                passing_variants.append((variant, flag, info))
                passing_variants.append((variant, flag, info))

        sys.stdout.write("Writing {} passing variants (of {}) to sample report\n".format(len(passing_variants),
                                                                                         iterated))
        utils.write_sample_variant_report(args.report, sample, passing_variants, args.variant_callers, thresholds)

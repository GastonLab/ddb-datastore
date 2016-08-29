#!/usr/bin/env python

import argparse
import getpass
import sys

from cassandra.auth import PlainTextAuthProvider
from cassandra.cqlengine import connection
from ddb import configuration

import utils
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

    thresholds = {'min_saf': 0.000001,
                  'max_maf': 0.005,
                  'depth': 200,
                  'regions': config['actionable_regions']}

    sys.stdout.write("Processing samples\n")
    for sample in samples:
        sys.stdout.write("Running Cassandra query for sample {}\n".format(sample))
        variants = SampleVariant.objects.timeout(None).filter(SampleVariant.reference_genome == config['genome_version'],
                                                              SampleVariant.sample == samples[sample]['sample_name'],
                                                              SampleVariant.run_id == samples[sample]['run_id'],
                                                              SampleVariant.library_name == samples[sample]['library_name'],
                                                              SampleVariant.max_som_aaf >= thresholds['min_saf'],
                                                              SampleVariant.max_maf_all <= thresholds['max_maf'],
                                                              SampleVariant.max_depth >= thresholds['depth']
                                                              ).allow_filtering()

        ordered_variants = variants.order_by('library_name', 'chr', 'pos',
                                             'ref', 'alt').limit(variants.count() + 1000)

        # variant_coords = samples[sample]['variant_coords'].split(',')

        validation_variants = list()
        for variant in ordered_variants:
            if 'cosmic_ids' in samples[sample].keys():
                cosmic_ids = samples[sample]['cosmic_ids'].split(',')
                if variant.cosmic_ids:
                    for cosmic_id in variant.cosmic_ids:
                        if cosmic_id in cosmic_ids:
                            validation_variants.append(variant)
                            break

            if 'rs_ids' in samples[sample].keys():
                rs_ids = samples[sample]['rs_ids'].split(',')
                if variant.rs_ids:
                    for rs_id in variant.rs_ids:
                        if rs_id in rs_ids:
                            validation_variants.append(variant)
                            break

            if 'amplicons' in samples[sample].keys():
                overlapping_amplicons = variant.amplicon_data['amplicon'].split(',')
                amplicons = samples[sample]['amplicons'].split(',')
                if variant.amplicon_data['amplicon']:
                    for amplicon in overlapping_amplicons:
                        if amplicons in amplicons:
                            validation_variants.append(variant)
                            break

        sys.stdout.write("Retrieved {} total variants\n".format(variants.count()))
        sys.stdout.write("Writing {} variants to sample report\n".format(len(validation_variants)))
        utils.write_sample_variant_report(args.report, sample, validation_variants, args.variant_callers)
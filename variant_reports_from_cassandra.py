#!/usr/bin/env python

import argparse
import getpass
import sys

from cassandra.auth import PlainTextAuthProvider
from cassandra.cqlengine import connection
from ddb import configuration

import utils
from variantstore import Variant

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

    thresholds = {'min_saf': 0.00000000001,
                  'max_maf': 0.005,
                  'depth': 200}

    sys.stdout.write("Processing samples\n")

    sys.stdout.write("Running Cassandra query\n")
    output_variants = list()
    for sample in samples:
        variants = Variant.objects.timeout(None).filter(Variant.reference_genome == config['genome_version'],
                                                        Variant.sample == samples[sample]['sample_name'],
                                                        Variant.run_id == samples[sample]['run_id'],
                                                        Variant.library_name == samples[sample]['library_name'],
                                                        ).allow_filtering()

        sys.stdout.write("Retrieved {} total variants\n".format(variants.count()))
        output_variants.extend(variants)

    utils.write_variant_report(args.report, output_variants, args.variant_callers)

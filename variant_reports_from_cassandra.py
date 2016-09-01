#!/usr/bin/env python

import argparse
import getpass
import sys
import csv

from cassandra.auth import PlainTextAuthProvider
from cassandra.cqlengine import connection
from ddb import configuration

import utils
from variantstore import Variant


def get_variants_list(infile):
    variants_list = list()

    with open(infile, 'r') as var_file:
        reader = csv.reader(var_file)
        for row in reader:
            var = dict()
            variants_list.append(var)

    return variants_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)")
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    parser.add_argument('-v', '--variant_callers', help="Comma-delimited list of variant callers used")
    parser.add_argument('-d', '--depth', help='Depth threshold', default=200)
    parser.add_argument('-m', '--max_pop_freq', help='Maximum population frequency threshold', default=0.005)
    parser.add_argument('-f', '--min_somatic_allele_freq', help='Minimum somatic frequency threshold', default=0.01)

    args = parser.parse_args()

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    variants = get_variants_list()

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "variantstore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "variantstore")

    thresholds = {'min_saf': args.min_somatic_allele_freq,
                  'max_maf': args.max_pop_freq,
                  'depth': args.depth}

    sys.stdout.write("Processing samples\n")

    sys.stdout.write("Running Cassandra queries\n")
    output_variants = list()
    for variant in variants:
        sample_variants = Variant.objects.timeout(None).filter(Variant.reference_genome == config['genome_version'],
                                                               Variant.chr == variant['chr'],
                                                               Variant.pos == variant['pos'],
                                                               Variant.ref == variant['ref']
                                                               ).allow_filtering()

    utils.write_variant_report(args.report, output_variants, args.variant_callers)

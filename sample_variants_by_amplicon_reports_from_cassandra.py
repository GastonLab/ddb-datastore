#!/usr/bin/env python

import argparse
import getpass
import sys
import csv

from cassandra.auth import PlainTextAuthProvider
from cassandra.cqlengine import connection
from ddb import configuration

import utils
from variantstore import TargetVariant
from collections import defaultdict


def get_amplicons_list(infile):
    amplicons_list = list()

    with open(infile, 'r') as amp_file:
        reader = csv.reader(amp_file)
        for row in reader:
            amplicons_list.append(row[0])

    return amplicons_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-l', '--list', help="Amplicon list file")
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

    amplicons = get_amplicons_list(args.list)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "variantstore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "variantstore")

    thresholds = {'min_saf': args.min_somatic_allele_freq,
                  'max_maf': args.max_pop_freq,
                  'depth': args.depth}
    callers = args.variant_callers.split(',')

    sys.stdout.write("Processing amplicons\n")

    sys.stdout.write("Running Cassandra queries\n")
    sample_variants = defaultdict(list)
    for amplicon in amplicons:
        target_variants = TargetVariant.objects.timeout(None).filter(TargetVariant.target == amplicon,
                                                                     TargetVariant.reference_genome == config['genome_version'],
                                                                     TargetVariant.max_som_aaf >= thresholds['min_saf'],
                                                                     TargetVariant.max_maf_all <= thresholds['max_maf'],
                                                                     TargetVariant.max_depth >= thresholds['depth']
                                                                     ).allow_filtering()

        ordered_variants = target_variants.order_by('sample', 'library_name', 'run_id', 'chr', 'pos',
                                                    'ref', 'alt').limit(target_variants.count() + 1000)
        for variant in ordered_variants:
            for caller in callers:
                if caller in variant.callers:
                    sample_variants[variant.sample].extend(ordered_variants)
                    break

    for sample in sample_variants:
        report_name = "{}.{}.txt".format(sample, args.report)
        utils.write_amplicon_variant_report(report_name, sample_variants[sample], args.variant_callers)

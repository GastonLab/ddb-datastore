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


def get_samples_list(infile):
    samples_list = list()

    with open(infile, 'r') as samp_file:
        reader = csv.reader(samp_file)
        for row in reader:
            samples_list.append(row[0])

    return samples_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--list', help="Amplicon list file")
    parser.add_argument('-s', '--samples', help="List of sample IDs")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)")
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    parser.add_argument('-v', '--variant_callers', help="Comma-delimited list of variant callers used")
    parser.add_argument('-d', '--depth', help='Depth threshold', default=200)
    parser.add_argument('-m', '--max_pop_freq', help='Maximum population frequency threshold', default=0.005)
    parser.add_argument('-f', '--min_somatic_allele_freq', help='Minimum somatic frequency threshold', default=0.01)

    args = parser.parse_args()

    amplicons = get_amplicons_list(args.list)
    samples = list()
    if args.samples:
        samples = get_samples_list(args.samples)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "variantstore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "variantstore")

    callers = args.variant_callers.split(',')

    sys.stdout.write("Processing amplicons\n")

    sys.stdout.write("Running Cassandra queries\n")
    for amplicon in amplicons:
        sys.stdout.write("Running query for amplicon: {}\n".format(amplicon))
        target_variants = TargetVariant.objects.timeout(None).filter(TargetVariant.target == amplicon,
                                                                     TargetVariant.reference_genome == 'GRCh37.75'
                                                                     ).allow_filtering()

        ordered_variants = target_variants.order_by('sample', 'library_name', 'run_id', 'chr', 'pos',
                                                    'ref', 'alt').limit(target_variants.count() + 1000)
        filtered_variants = defaultdict(lambda: defaultdict(int))

        for variant in ordered_variants:
            for caller in callers:
                if caller in variant.callers:
                    if args.samples:
                        if variant.sample in samples:
                            filtered_variants.append(variant)
                            break
                    else:
                        filtered_variants.append(variant)
                        break

        report_name = "{}.{}.txt".format(amplicon, args.report)
        utils.write_amplicon_variant_report(report_name, filtered_variants, args.variant_callers)

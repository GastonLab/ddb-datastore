#!/usr/bin/env python

import argparse
import getpass
import sys
import csv

from cassandra.auth import PlainTextAuthProvider
from cassandra.cqlengine import connection
from ddb import configuration

import utils
from coveragestore import AmpliconCoverage
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
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-l', '--list', help="Amplicon list file")
    parser.add_argument('-s', '--samples', help="List of sample IDs")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)")
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)

    args = parser.parse_args()

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

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

    sys.stdout.write("Processing amplicons\n")

    sys.stdout.write("Running Cassandra queries\n")
    sample_amplicons = defaultdict(list)
    for amplicon_name in amplicons:
        sys.stdout.write("Running queries for amplicon : {}\n".format(amplicon_name))
        target_coverage = AmpliconCoverage.objects.timeout(None).filter(amplicon=amplicon_name).allow_filtering()
        sys.stdout.write("Returned {} results\n".format(target_coverage.count()))

        ordered_coverage = target_coverage.order_by('sample').limit(target_coverage.count() + 10000)

        for amplicon in target_coverage:
            sys.stdout.write("Amplicon: {}\tSample:{}\n".format(amplicon.amplicon, amplicon.sample))
            if args.samples:
                if amplicon.sample in samples:
                    sample_amplicons[amplicon.sample].append(amplicon)
                    break
            else:
                sample_amplicons[amplicon.sample].append(amplicon)
                break

    for sample in sample_amplicons:
        report_name = "{}.{}.txt".format(sample, args.report)
        utils.write_amplicon_coverage_report(report_name, sample_amplicons[sample])

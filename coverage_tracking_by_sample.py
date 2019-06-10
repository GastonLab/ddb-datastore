#!/usr/bin/env python

import sys
import utils
import getpass
import argparse

import numpy as np
from collections import defaultdict
from coveragestore import AmpliconCoverage
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--list', help="File containing list of amplicon names to check")
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    args = parser.parse_args()
    args.logLevel = "INFO"

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "coveragestore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "coveragestore")

    amplicons_list = list()
    amplicon_coverage_stats = defaultdict(dict)
    amplicon_stats_by_month = defaultdict(lambda: defaultdict(dict))
    amplicon_stats_by_run = defaultdict(lambda: defaultdict(dict))

    sys.stdout.write("Parsing Amplicon List\n")
    target_amplicons = utils.get_target_amplicons(args.list)
    for amplicon in target_amplicons:
        if amplicon not in amplicons_list:
            amplicons_list.append(amplicon)

    sys.stdout.write("Processing Amplicon Data\n")
    for amplicon in target_amplicons:
        sys.stdout.write("Retrieving Coverage Data for {}\n".format(amplicon))
        coverage_values = list()
        coverage_data = AmpliconCoverage.objects.timeout(None).filter(
            AmpliconCoverage.amplicon == amplicon
        )

        ordered_samples = coverage_data.order_by('sample', 'run_id').limit(coverage_data.count() + 1000)
        sys.stderr.write("There are {} samples retrieved\n".format(coverage_data.count()))
        sys.stdout.write("Sample\tLibrary\tRunID\tCov\n")
        num = 0
        for result in ordered_samples:
            num += 1
            coverage_values.append(result.mean_coverage)
            sys.stdout.write("{}\t{}\t{}\t{}\n".format(result.sample,
                                                       result.library_name,
                                                       result.run_id,
                                                       result.mean_coverage))
        sys.stderr.write("Iterated {} times\n".format(num))

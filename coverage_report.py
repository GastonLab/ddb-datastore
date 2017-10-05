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
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)", default='report')
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

    target_amplicons = utils.get_target_amplicons(args.list)
    for amplicon in target_amplicons:
        if amplicon not in amplicons_list:
            amplicons_list.append(amplicon)

    for amplicon in target_amplicons:
        coverage_values = list()
        coverage_values_by_month = defaultdict(list)
        coverage_data = AmpliconCoverage.objects.timeout(None).filter(
            AmpliconCoverage.amplicon == amplicon
        )
        ordered_samples = coverage_data.order_by('sample', 'run_id').limit(coverage_data.count() + 1000)
        for result in ordered_samples:
            coverage_values.append(result.mean_coverage)
            yr_month_id = result.run_id[:4]
            coverage_values_by_month[yr_month_id].append(result.mean_coverage)

        amplicon_coverage_stats[amplicon]['median'] = np.median(coverage_values)
        amplicon_coverage_stats[amplicon]['std_dev'] = np.std(coverage_values)
        amplicon_coverage_stats[amplicon]['min'] = np.amin(coverage_values)
        amplicon_coverage_stats[amplicon]['max'] = np.amax(coverage_values)

        for yr_month_id in coverage_values_by_month:
            amplicon_stats_by_month[amplicon][yr_month_id]['median'] = np.median(coverage_values)
            amplicon_stats_by_month[amplicon][yr_month_id]['std_dev'] = np.std(coverage_values)
            amplicon_stats_by_month[amplicon][yr_month_id]['min'] = np.amin(coverage_values)
            amplicon_stats_by_month[amplicon][yr_month_id]['max'] = np.amax(coverage_values)

    with open("coverage_analysis.txt", "w") as coverage_report:
        coverage_report.write("Amplicon\tMedian\Std\tmin\tmax\n")
        for amplicon in amplicon_coverage_stats:
            coverage_report.write("{}\t{}\t{}\t{}\t{}\n".format(amplicon, amplicon_coverage_stats[amplicon]['median'],
                                                                amplicon_coverage_stats[amplicon]['std_dev'],
                                                                amplicon_coverage_stats[amplicon]['min'],
                                                                amplicon_coverage_stats[amplicon]['max']))

    with open("coverage_analysis_by_month.txt", "w") as coverage_report:
        coverage_report.write("Amplicon\tMonth\tMedian\Std\tmin\tmax\n")
        for amplicon in amplicon_stats_by_month:
            for yr_month_id in amplicon_stats_by_month[amplicon]:
                coverage_report.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(amplicon, yr_month_id,
                                                                        amplicon_stats_by_month[amplicon][yr_month_id]['median'],
                                                                        amplicon_stats_by_month[amplicon][yr_month_id]['std_dev'],
                                                                        amplicon_stats_by_month[amplicon][yr_month_id]['min'],
                                                                        amplicon_stats_by_month[amplicon][yr_month_id]['max']))

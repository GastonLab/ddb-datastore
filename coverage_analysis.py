#!/usr/bin/env python

import sys
import csv
import argparse
import utils
import getpass
import plotly
import plotly.plotly as py
import plotly.graph_objs as go


from collections import defaultdict

from coveragestore import AmpliconCoverage
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider


def get_regions(infile):
    regions = list()
    with open(infile, 'r') as csvfile:
        reader = csv.reader(csvfile, dialect='excel-tab')
        for row in reader:
            regions.append(row[3])

    return regions


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--regions', help="Regions file with Amplicon names")
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-n', '--num_libraries', help='Number of libraries to consider when calculating coverage '
                                                      'stats', default=8)
    parser.add_argument('-s', '--stat', help="Statistic to plot")
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)

    args = parser.parse_args()

    regions = get_regions(args.regions)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "coveragestore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "coveragestore")

    summary_data = defaultdict(lambda: defaultdict(list))
    sys.stdout.write("Calculation coverage data per amplicon/region\n")
    num_regions = len(regions)
    for region in regions:
        sys.stdout.write("Running Cassandra query for region {}\n".format(region))
        samples = AmpliconCoverage.objects.timeout(None).filter(amplicon=region,
                                                                num_libraries_in_run=args.num_libraries).allow_filtering()

        ordered_samples = samples.order_by('sample', 'run_id', 'library_name').limit(samples.count() + 1000)

        # sys.stdout.write("Retrieved data for {} libraries\n".format(samples.count()))
        passing_variants = list()
        passed = 0
        iterated = 0
        data = defaultdict(lambda: defaultdict(list))
        counts = defaultdict(int)
        for sample in ordered_samples:
            counts[sample.extraction] += 1
            data[sample.extraction]['num_reads'].append(sample.num_reads)
            data[sample.extraction]['mean_coverage'].append(sample.mean_coverage)
            data[sample.extraction]['percent_bp1'].append(sample.perc_bp_cov_at_thresholds[500])
            data[sample.extraction]['percent_bp2'].append(sample.perc_bp_cov_at_thresholds[1000])

            summary_data[sample.extraction]['num_reads'].append(sample.num_reads)
            summary_data[sample.extraction]['mean_coverage'].append(sample.mean_coverage)
            summary_data[sample.extraction]['percent_bp1'].append(sample.perc_bp_cov_at_thresholds[500])
            summary_data[sample.extraction]['percent_bp2'].append(sample.perc_bp_cov_at_thresholds[1000])

            iterated += 1

        # sys.stdout.write("Processed {} libraries\n".format(iterated))
        # traces = list()
        # for extraction in data:
        #     # sys.stdout.write("Extraction: {}, Num Libraries {}\n".format(extraction, counts[extraction]))
        #     trace = go.Box(y=data[extraction][args.stat], boxpoints='all', jitter=0.3, pointpos=1.8,
        #                    name="{} ({})".format(extraction, len(data[extraction][args.stat])))
        #     traces.append(trace)
        #
        # plotly.offline.plot(traces, filename="{}_{}_boxplot.html".format(region, args.stat))

    traces = list()
    for extraction in summary_data:
        # sys.stdout.write("Extraction: {}, Num Libraries {}\n".format(extraction, counts[extraction]))
        trace = go.Box(y=summary_data[extraction][args.stat], boxpoints='all', jitter=0.3, pointpos=1.8,
                       name="{} ({})".format(extraction, len(summary_data[extraction][args.stat]) / num_regions))
        traces.append(trace)

    plotly.offline.plot(traces, filename="{}_{}_{}_regions-boxplot.html".format("All", args.stat, num_regions))

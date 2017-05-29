#!/usr/bin/env python

import argparse
import csv
import getpass
import sys
from collections import defaultdict

import plotly
import plotly.graph_objs as go
from cassandra.auth import PlainTextAuthProvider
from cassandra.cqlengine import connection

from ddb_data.coveragestore import AmpliconCoverage


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
    parser.add_argument('-i', '--include_simulated', help="Statistic to plot", default=False)
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    parser.add_argument('-o', '--output', help="Output file name for summary of problematic data",
                        default="coverage_issues.txt",)

    args = parser.parse_args()

    regions = get_regions(args.regions)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "coveragestore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "coveragestore")

    summary_data = defaultdict(lambda: defaultdict(list))
    problem_amplicons1 = defaultdict(list)
    problem_amplicons2 = defaultdict(list)
    problem_amplicon_traces = list()
    sys.stdout.write("Calculation coverage data per amplicon/region\n")
    num_regions = len(regions)

    with open(args.output, 'w') as output:
        output.write("Amplicon\tNum Samples\tNum Samples < 10% bp @ 500x\tNum Samples < 10% bp @ 1000x\t"
                     "Num Samples < 20% bp @ 500x\tNum Samples < 20% bp @ 1000x\tNum Samples < 30% bp @ 500x\t"
                     "Num Samples < 30% bp @ 1000x\tNum Samples < 50% bp @ 500x\tNum Samples < 50% bp @ 1000x\t"
                     "Num Samples < 80% bp @ 500x\tNum Samples < 80% bp @ 1000x\tNum Samples < 100% bp @ 500x\t"
                     "Num Samples < 100% bp @ 1000x\tNum Samples OK 500x\tNum Samples OK 1000x\n")
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
            problem_counts = defaultdict(lambda: defaultdict(int))
            for sample in ordered_samples:
                if not args.include_simulated:
                    if sample.library_name.startswith("subsample-"):
                        continue

                counts[sample.extraction] += 1
                data[sample.extraction]['num_reads'].append(sample.num_reads)
                data[sample.extraction]['mean_coverage'].append(sample.mean_coverage)
                data[sample.extraction]['percent_bp1'].append(sample.perc_bp_cov_at_thresholds[500])
                data[sample.extraction]['percent_bp2'].append(sample.perc_bp_cov_at_thresholds[1000])

                summary_data[sample.extraction]['num_reads'].append(sample.num_reads)
                summary_data[sample.extraction]['mean_coverage'].append(sample.mean_coverage)
                summary_data[sample.extraction]['percent_bp1'].append(sample.perc_bp_cov_at_thresholds[500])
                summary_data[sample.extraction]['percent_bp2'].append(sample.perc_bp_cov_at_thresholds[1000])

                if sample.perc_bp_cov_at_thresholds[500] < 100.00:
                    problem_amplicons1[region].append((sample.sample, sample.run_id, sample.library_name))

                if sample.perc_bp_cov_at_thresholds[1000] < 100.00:
                    problem_amplicons2[region].append((sample.sample, sample.run_id, sample.library_name))

                if sample.perc_bp_cov_at_thresholds[500] <= 10.00:
                    problem_counts[500][10] += 1
                elif 10.00 < sample.perc_bp_cov_at_thresholds[500] <= 20.00:
                    problem_counts[500][20] += 1
                elif 20.00 < sample.perc_bp_cov_at_thresholds[500] <= 30.00:
                    problem_counts[500][30] += 1
                elif 30.00 < sample.perc_bp_cov_at_thresholds[500] <= 50.00:
                    problem_counts[500][50] += 1
                elif 50.00 < sample.perc_bp_cov_at_thresholds[500] <= 80.00:
                    problem_counts[500][80] += 1
                elif 80.00 < sample.perc_bp_cov_at_thresholds[500] < 100.00:
                    problem_counts[500][100] += 1
                else:
                    problem_counts[500]['ok'] += 1

                if sample.perc_bp_cov_at_thresholds[1000] <= 10.00:
                    problem_counts[1000][10] += 1
                elif 10.00 < sample.perc_bp_cov_at_thresholds[1000] <= 20.00:
                    problem_counts[1000][20] += 1
                elif 20.00 < sample.perc_bp_cov_at_thresholds[1000] <= 30.00:
                    problem_counts[1000][30] += 1
                elif 30.00 < sample.perc_bp_cov_at_thresholds[1000] <= 50.00:
                    problem_counts[1000][50] += 1
                elif 50.00 < sample.perc_bp_cov_at_thresholds[1000] <= 80.00:
                    problem_counts[1000][80] += 1
                elif 80.00 < sample.perc_bp_cov_at_thresholds[1000] < 100.00:
                    problem_counts[1000][100] += 1
                else:
                    problem_counts[1000]['ok'] += 1

                iterated += 1

            if problem_counts[500]['ok'] < iterated:
                output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
                             "\n".format(region, iterated,
                                         problem_counts[500][10], problem_counts[1000][10],
                                         problem_counts[500][20], problem_counts[1000][20],
                                         problem_counts[500][30], problem_counts[1000][30],
                                         problem_counts[500][50], problem_counts[1000][50],
                                         problem_counts[500][80], problem_counts[1000][80],
                                         problem_counts[500][100], problem_counts[1000][100],
                                         problem_counts[500]['ok'], problem_counts[1000]['ok']))

                for extraction in data:
                    trace = go.Box(y=data[extraction][args.stat], boxpoints='all', jitter=0.3, pointpos=1.8,
                                   name="{} ({})".format(region, extraction))
                    problem_amplicon_traces.append(trace)

    boxplots = list()
    histograms = list()
    for extraction in summary_data:
        # sys.stdout.write("Extraction: {}, Num Libraries {}\n".format(extraction, counts[extraction]))
        boxplot = go.Box(y=summary_data[extraction][args.stat], boxpoints='all', jitter=0.3, pointpos=1.8,
                         name="{} ({})".format(extraction, len(summary_data[extraction][args.stat]) / num_regions))
        histogram = go.Histogram(x=summary_data[extraction][args.stat], opacity=0.75, autobinx=False,
                                 xbins=dict(
                                    start=0,
                                    end=100,
                                    size=10
                                 ))

        boxplots.append(boxplot)
        histograms.append(histogram)

    plotly.offline.plot(boxplots, filename="{}_{}_{}_regions-boxplot.html".format("All", args.stat, num_regions))

    layout = go.Layout(
        barmode='overlay'
    )

    fig = go.Figure(data=histograms, layout=layout)
    plotly.offline.plot(fig, filename="{}_{}_{}_regions-distribution.html".format("All", args.stat, num_regions))

    plotly.offline.plot(problem_amplicon_traces, filename="problem_amplicons_{}.html".format(args.stat))

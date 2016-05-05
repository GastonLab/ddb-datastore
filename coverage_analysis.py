#!/usr/bin/env python

import sys
import csv
import argparse
import utils
import getpass
import plotly

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
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)

    args = parser.parse_args()

    regions = get_regions(args.regions)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "variantstore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "variantstore")

    thresholds = {'threshold1': 500,
                  'threshold2': 1000}

    sys.stdout.write("Calculation coverage data per amplicon/region\n")
    for region in regions:
        sys.stdout.write("Running Cassandra query for region {}\n".format(region))
        samples = AmpliconCoverage.objects.timeout(None).filter(amplicon=region,
                                                                num_libraries_in_run=args.num_libraries).allow_filtering()

        ordered_samples = samples.order_by('sample', 'run_id', 'library_name').limit(samples.count() + 1000)

        sys.stdout.write("Retrieved data for {} libraries\n".format(samples.count()))
        passing_variants = list()
        passed = 0
        iterated = 0
        data = defaultdict(list)
        counts = defaultdict(int)
        for sample in ordered_samples:
            counts[sample.extraction] += 1
            data[sample.extraction].append(sample.num_reads)

            iterated += 1

        sys.stdout.write("Processed {} libraries\n".format(iterated))
        for extraction in data:
            sys.stdout.write("Extraction: {}, Num Libraries {}\n".format(extraction, counts[extraction]))

        # plotly.offline.plot()

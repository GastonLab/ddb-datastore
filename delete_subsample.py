#!/usr/bin/env python

import sys
import csv
import argparse
import getpass

from collections import defaultdict

from coveragestore import AmpliconCoverage

from cassandra.cluster import Cluster
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
        connection.setup([args.address], "coveragestore", auth_provider=auth_provider)
        cluster = Cluster([args.address], auth_provider=auth_provider)
    else:
        connection.setup([args.address], "coveragestore")
        cluster = Cluster([args.address])

    session = cluster.connect(keyspace="coveragestore")

    summary_data = defaultdict(lambda: defaultdict(list))
    problem_amplicons1 = defaultdict(list)
    problem_amplicons2 = defaultdict(list)
    problem_amplicon_traces = list()
    sys.stdout.write("Calculation coverage data per amplicon/region\n")
    num_regions = len(regions)

    for region in regions:
        sys.stdout.write("Running Cassandra query for region {}\n".format(region))
        samples = AmpliconCoverage.objects.timeout(None).filter(amplicon=region,
                                                                num_libraries_in_run=args.num_libraries).allow_filtering()

        ordered_samples = samples.order_by('sample', 'run_id', 'library_name').limit(samples.count() + 1000)

        passing_variants = list()
        passed = 0
        iterated = 0
        data = defaultdict(lambda: defaultdict(list))
        counts = defaultdict(int)
        problem_counts = defaultdict(lambda: defaultdict(int))
        for sample in ordered_samples:
            if sample.library_name.startswith("subsample-"):
                counts[sample.extraction] += 1
                session.execute("DELETE FROM amplicon_coverage WHERE amplicon = '{}' AND sample = '{}' AND "
                                "run_id = '{}' AND library_name = '{}' AND program_name = '{}'"
                                "".format(region, sample.sample, sample.run_id, sample.library_name,
                                          sample.program_name))
                session.execute("DELETE FROM sample_coverage WHERE sample = '{}' AND amplicon = '{}' AND "
                                "run_id = '{}' AND library_name = '{}' AND program_name = '{}'"
                                "".format(region, sample.sample, sample.run_id, sample.library_name,
                                          sample.program_name))
                iterated += 1

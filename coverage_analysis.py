#!/usr/bin/env python

import sys
import csv
import argparse
import utils
import getpass

from coveragestore import AmpliconCoverage
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider


def get_regions(infile):
    regions = list()
    with open(infile, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            regions.append(row[3])

    return regions


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--regions', help="Regions file with Amplicon names")
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
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
        samples = AmpliconCoverage.objects.timeout(None).filter(amplicon=region).allow_filtering()

        sys.stdout.write("Retrieved data for {} total samples\n".format(samples.count()))
        sys.stdout.write("Running filters on sample variants\n")
        passing_variants = list()
        passed = 0
        iterated = 0
        for sample in samples:
            iterated += 1

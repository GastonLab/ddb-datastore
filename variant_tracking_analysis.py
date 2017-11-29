#!/usr/bin/env python

import sys
import csv
import utils
import getpass
import argparse

import numpy as np
from collections import defaultdict
from variantstore import Variant
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--list', help="File containing list of variants to check")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)", default='report')
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    args = parser.parse_args()
    args.logLevel = "INFO"

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "variantstore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "variantstore")

    sys.stdout.write("Proccessng through selected variants")
    with open(args.list, "r") as variants_list:
        with open("variant_analysis_tracking_{}.txt".format(args.report), "w") as output:
            output.write("Chr\tPos\tRef\tAlt\tCodon\tAA\tNum in DB\tMedian VAF\tStd VAF\tAmplicon\n")
            reader = csv.reader(variants_list, dialect='excel-tab')
            for row in reader:
                match_variants = Variant.objects.timeout(None).filter(
                    Variant.reference_genome == "GRCh37.75",
                    Variant.chr == row[0],
                    Variant.pos == row[1],
                    Variant.ref == row[3],
                    Variant.alt == row[4]
                ).allow_filtering()

                num_matches = match_variants.count()
                ordered_var = match_variants.order_by('ref', 'alt', 'sample', 'library_name',
                                                      'run_id').limit(num_matches + 1000)
                vafs = list()
                run_vafs = list()
                num_times_callers = defaultdict(int)
                seen_libraries = list()

                variant_data = dict()

                for var in ordered_var:
                    if var.library_name not in seen_libraries:
                        seen_libraries.append(var.library_name)
                        vaf = var.max_som_aaf
                        vafs.append(vaf)
                        variant_data['codon'] = var.codon_change
                        for caller in var.callers:
                            num_times_callers[caller] += 1

                vaf_median = np.median(vafs)
                vaf_std_dev = np.std(vafs)

                caller_counts_elements = list()
                for caller in num_times_callers:
                    caller_counts_elements.append("{}: {}".format(caller, num_times_callers[caller]))
                num_times_callers = ",".join(caller_counts_elements)
                output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t"
                             "\n".format(row[0], row[1], row[3], row[4], row[5], row[6], num_matches, vaf_median,
                                         vaf_std_dev, row[7]))


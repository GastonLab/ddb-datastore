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
    parser.add_argument('-l', '--list', help="File containing list of variants\
                        to check")
    parser.add_argument('-r', '--report', help="Root name for reports\
                        (per sample)", default='report')
    parser.add_argument('-a', '--address', help="IP Address for Cassandra\
                        connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for\
                        login', default=None)
    args = parser.parse_args()
    args.logLevel = "INFO"

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username,
                                              password=password)
        connection.setup([args.address], "variantstore",
                         auth_provider=auth_provider)
    else:
        connection.setup([args.address], "variantstore")

    sys.stdout.write("Proccessng through selected variants")
    with open(args.list, "r") as variants_list:
        with open("variant_pos_samples_{}.txt".format(args.report), "w") as output:
            output.write("Chr\tPos\tRef\tAlt\tCodon\tAA\tAmplicon\tSample\tLibrary\tRun\tVAF\tCallers\n")
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

                for var in ordered_var:
                    output.write("{chr}\t{pos}\t{ref}\t{alt}\t{codon}\t{aa}\t{amplicon}\t{sample}\t{lib}\t{run}\t{vaf}\t{call}"
                                 "\n".format(chr=var.chr, pos=var.pos, ref=var.ref, alt=var.alt, codon=var.codon_change,
                                             aa=var.aa_change, amplicon=var.amplicon, sample=var.sample, lib=var.library_name,
                                             run=var.run_id, vaf=var.max_som_aaf, call=",".join(var.callers) or None))

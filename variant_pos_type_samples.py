#!/usr/bin/env python

import os
import sys
import csv
import fnmatch
import getpass
import argparse

from ddb import configuration
from collections import defaultdict
from variantstore import Variant
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--list', help="File containing list of variants\
                        to check")
    parser.add_argument('-r', '--report', help="Root name for report\
                        (per sample)", default='report')
    parser.add_argument('-a', '--address', help="IP Address for Cassandra\
                        connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for\
                        login', default=None)
    parser.add_argument('-c', '--configuration',
                        help="Configuration file for various settings")
    args = parser.parse_args()
    args.logLevel = "INFO"

    type = "colorectal"

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username,
                                              password=password)
        connection.setup([args.address], "variantstore",
                         auth_provider=auth_provider)
    else:
        connection.setup([args.address], "variantstore")

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    type_samples = list()

    for root, dirs, files in os.walk("."):
        for samples_file in fnmatch.filter(files, "1*_M0373?.config"):
            sys.stderr.write("Reading file: {}\n".format(os.path.join(root, samples_file)))

            sys.stdout.write("Parsing sample data\n")
            libraries = configuration.configure_samples(os.path.join(root, samples_file), config)
            samples = configuration.merge_library_configs_samples(libraries)
            for sample in samples:
                for library in samples[sample]:
                    if samples[sample][library]['report'].startswith(type):
                        print "Colorectal case sequencing library found: {}\n".format(samples[sample][library]['library_name'])
                        type_samples.append(samples[sample][library]['library_name'])

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

                for var in ordered_var:
                    if(var.library_name in type_samples):
                        output.write("{chr}\t{pos}\t{ref}\t{alt}\t{codon}\t{aa}\t{amplicon}\t{sample}\t{lib}\t{run}\t{vaf}\t{call}"
                                     "\n".format(chr=var.chr, pos=var.pos, ref=var.ref, alt=var.alt, codon=var.codon_change,
                                                 aa=var.aa_change, amplicon=var.amplicon, sample=var.sample, lib=var.library_name,
                                                 run=var.run_id, vaf=var.max_som_aaf, call=",".join(var.callers) or None))

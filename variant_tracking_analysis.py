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
    parser.add_argument('-r', '--report', help="Root name for reports", default='report')
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

    sys.stdout.write("Proccessng through selected variants\n")
    with open(args.list, "r") as variants_list:
        with open("variant_analysis_tracking_{}.txt".format(args.report), "w") as output:
            output.write("Sample\tLibrary\tGene\tAmplicon\tRef\tAlt\tCodon\tAA\t"
                         "max_somatic_aaf\tCallers\tCOSMIC_IDs\tCOSMIC_NumSamples\tCOSMIC_AA\t"
                         "Clin_Sig\tClin_HGVS\tClin_Disease\t"
                         "Coverage\tNum Reads\tImpact\tSeverity\tmax_maf_all\tmax_maf_no_fin\t"
                         "min_caller_depth\tmax_caller_depth\tChrom\tStart\tEnd\trsIDs\n")
            reader = csv.reader(variants_list, dialect='excel-tab')
            for row in reader:
                # sys.stdout.write("Processing row: {}\n".format(row))
                match_variants = Variant.objects.timeout(None).filter(
                    Variant.reference_genome == "GRCh37.75",
                    Variant.chr == row[0],
                    Variant.pos == row[1],
                    Variant.ref == row[2],
                    Variant.alt == row[3]
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
                    output.write("{sample}\t{library}\t{gene}\t{amp}\t{ref}\t{alt}\t{codon}\t{aa}\t"
                                 "{max_som_aaf}\t{callers}\t{cosmic}\t{cosmic_nsamples}\t{cosmic_aa}\t{csig}\t"
                                 "{hgvs}\t{cdis}\t{impact}\t{severity}\t{max_maf_all}\t{max_maf_no_fin}\t"
                                 "{min_depth}\t{max_depth}\t{chr}\t{start}\t{end}\t{rsids}"
                                 "\n".format(sample=var.sample,
                                               library=var.library_name,
                                               chr=var.chr,
                                               start=var.pos,
                                               end=var.end,
                                               gene=var.gene,
                                               ref=var.ref,
                                               alt=var.alt,
                                               codon=var.codon_change,
                                               aa=var.aa_change,
                                               rsids=",".join(var.rs_ids),
                                               cosmic=",".join(var.cosmic_ids) or None,
                                               cosmic_nsamples=var.cosmic_data['num_samples'],
                                               cosmic_aa=var.cosmic_data['aa'],
                                               amp=var.amplicon_data['amplicon'],
                                               csig=var.clinvar_data['significance'],
                                               hgvs=var.clinvar_data['hgvs'],
                                               cdis=var.clinvar_data['disease'],
                                               impact=var.impact,
                                               severity=var.severity,
                                               max_maf_all=var.max_maf_all,
                                               max_maf_no_fin=var.max_maf_no_fin,
                                               max_som_aaf=var.max_som_aaf,
                                               min_depth=var.min_depth,
                                               max_depth=var.max_depth,
                                               callers=",".join(var.callers) or None))

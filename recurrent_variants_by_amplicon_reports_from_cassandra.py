#!/usr/bin/env python

import argparse
import getpass
import sys
import csv

from cassandra.auth import PlainTextAuthProvider
from cassandra.cqlengine import connection
from ddb import configuration

import utils
from variantstore import TargetVariant
from collections import defaultdict


def get_amplicons_list(infile):
    amplicons_list = list()

    with open(infile, 'r') as amp_file:
        reader = csv.reader(amp_file)
        for row in reader:
            amplicons_list.append(row[0])

    return amplicons_list


def get_samples_list(infile):
    samples_list = list()

    with open(infile, 'r') as samp_file:
        reader = csv.reader(samp_file)
        for row in reader:
            samples_list.append(row[0])

    return samples_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--list', help="Amplicon list file")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)")
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    parser.add_argument('-m', '--max_pop_freq', help='Maximum population frequency threshold', default=0.005)
    parser.add_argument('-o', '--output', help="Output file name for variants")

    args = parser.parse_args()

    amplicons = get_amplicons_list(args.list)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "variantstore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "variantstore")

    sys.stdout.write("Running Cassandra queries\n")
    for amplicon in amplicons:
        sys.stdout.write("Running query for amplicon: {}\n".format(amplicon))
        target_variants = TargetVariant.objects.timeout(None).filter(TargetVariant.target == amplicon,
                                                                     TargetVariant.reference_genome == 'GRCh37.75'
                                                                     ).allow_filtering()

        ordered_variants = target_variants.order_by('sample', 'library_name', 'run_id', 'chr', 'pos',
                                                    'ref', 'alt').limit(target_variants.count() + 1000)

        filtered_variants = defaultdict(lambda: defaultdict(int))

        for variant in ordered_variants:
            if variant.max_maf_all <= args.max_pop_freq:
                key = "{}-{}-{}-{}-{}".format(variant.reference_genome, variant.chr, variant.pos, variant.ref,
                                              variant.alt)

                filtered_variants[key]['instances'] += 1
                filtered_variants[key]['panel'] = variant.panel_name
                filtered_variants[key]['amplicon'] = variant.amplicon_data['amplicon']
                filtered_variants[key]['gene'] = variant.gene or None
                filtered_variants[key]['amino_acid'] = variant.aa_change or None
                filtered_variants[key]['somatic_freq'] = variant.min_som_aaf
                filtered_variants[key]['cosmic_ids'] = ",".join(variant.cosmic_ids) or None
                filtered_variants[key]['num_cosmic'] = variant.cosmic_data['num_samples']
                # filtered_variants[key]['depth'] += variant.depth this doesn't exist yet

                if 'freebayes' in variant.callers:
                    filtered_variants[key]['freebayes'] += 1
                if 'scalpel' in variant.callers:
                    filtered_variants[key]['scalpel'] += 1
                if 'mutect' in variant.callers:
                    filtered_variants[key]['mutect'] += 1
                if 'pindel' in variant.callers:
                    filtered_variants[key]['pindel'] += 1
                if 'vardict' in variant.callers:
                    filtered_variants[key]['vardict'] += 1
                if 'platypus' in variant.callers:
                    filtered_variants[key]['platypus'] += 1

        with open(args.output, 'a') as outfile:
            outfile.write("Key\tGene\tAmplicon\tAmino Acid\tNum Instances\tNum Cosmic\tCosmic IDs\tSomatic AF\t"
                          "FreeBayes\tScalpel\tMuTect\tPindel\tVarDict\tPlatypus\n")
            for variant in filtered_variants.keys():
                outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
                              "\n".format(variant, filtered_variants[variant]['gene'],
                                          filtered_variants[variant]['amplicon'],
                                          filtered_variants[variant]['amino_acid'],
                                          filtered_variants[variant]['num_cosmic'],
                                          filtered_variants[variant]['cosmic_ids'],
                                          filtered_variants[variant]['somatic_freq'],
                                          filtered_variants[variant]['freebayes'],
                                          filtered_variants[variant]['scalpel'],
                                          filtered_variants[variant]['mutect'],
                                          filtered_variants[variant]['pindel'],
                                          filtered_variants[variant]['vardict'],
                                          filtered_variants[variant]['platypus']))

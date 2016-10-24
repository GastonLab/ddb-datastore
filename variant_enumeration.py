#!/usr/bin/env python

import argparse
import getpass
import sys

from cassandra.auth import PlainTextAuthProvider
from cassandra.cqlengine import connection

from variantstore import Variant
from collections import defaultdict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    parser.add_argument('-o', '--output', help='Output file name')

    args = parser.parse_args()

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "variantstore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "variantstore")

    variant_details = defaultdict(lambda: defaultdict(int))
    all_variants = Variant.objects.all()
    sys.stdout.write("Retrieved {} variants from the database\n".format(all_variants.count()))
    count = 0
    for variant in all_variants:
        count += 1
        key = "{}-{}-{}-{}-{}".format(variant.reference_genome, variant.chr, variant.pos, variant.ref, variant.alt)

        variant_details[key]['instances'] += 1
        variant_details[key]['panel'] = variant.panel_name
        variant_details[key]['gene'] = variant.gene
        variant_details[key]['amino_acid'] = variant.aa_change or None

        if 'freebayes' in variant.callers:
            variant_details[key]['freebayes'] += 1
        if 'scalpel' in variant.callers:
            variant_details[key]['scalpel'] += 1
        if 'mutect' in variant.callers:
            variant_details[key]['mutect'] += 1
        if 'pindel' in variant.callers:
            variant_details[key]['pindel'] += 1
        if 'vardict' in variant.callers:
            variant_details[key]['vardict'] += 1
        if 'platypus' in variant.callers:
            variant_details[key]['platypus'] += 1

    with open(args.output, 'w') as outfile:
        outfile.write("Variant Key\tGene\tAA\tPanel\tNum Instances\tNum MuTect\tNum FreeBayes\tNum VarDict\t"
                      "Num Platypus\tNum Scalpel\tNum Pindel\n")
        for variant in variant_details.keys():
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(variant,
                                                                                variant_details[variant]['gene'],
                                                                                variant_details[variant]['amino_acid'],
                                                                                variant_details[variant]['panel'],
                                                                                variant_details[variant]['instances'],
                                                                                variant_details[variant]['mutect'],
                                                                                variant_details[variant]['freebayes'],
                                                                                variant_details[variant]['vardict'],
                                                                                variant_details[variant]['platypus'],
                                                                                variant_details[variant]['scalpel'],
                                                                                variant_details[variant]['pindel']))

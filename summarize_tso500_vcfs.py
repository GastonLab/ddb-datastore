#!/usr/bin/env python

# This script is intended to act as a backup reporting system to convert VCF
# files directly to an excel report format, bypassing the VariantStore
# Cassandra database.

import re
import sys
import csv
import xlwt
import utils
import cyvcf2
import argparse

from cyvcf2 import VCF
from ddb import vcf_parsing
from ddb import configuration
from collections import defaultdict

def get_clinvar_info(variant):
    clinvar_data = dict()

    clinvar_data['significance'] = variant.INFO.get('clinvar_significance') or 'None'
    clinvar_data['pathogenic'] = variant.INFO.get('clinvar_pathogenic') or 'None'
    clinvar_data['hgvs'] = variant.INFO.get('clinvar_hgvs') or 'None'
    clinvar_data['revstatus'] = variant.INFO.get('clinvar_revstatus') or 'None'
    clinvar_data['org'] = variant.INFO.get('clinvar_org') or 'None'
    clinvar_data['disease'] = variant.INFO.get('clinvar_diseasename') or 'None'
    clinvar_data['accession'] = variant.INFO.get('clinvar_accession') or 'None'

    try:
        clinvar_data['origin'] = variant.INFO.get('clinvar_origin') or 'None'
    except IndexError:
        clinvar_data['origin'] = 'None'

    return clinvar_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',
                        help="Input configuration file for samples")
    parser.add_argument('-o', '--output',
                        help='Output file name for filtered variants',
                        default=0.005)
    args = parser.parse_args()
    args.logLevel = "INFO"

    with open(args.input, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        with open(args.output, 'w') as output:
            output.write("Sample\t")
            output.write("Gene\t")
            output.write("Ref\t")
            output.write("Alt\t")
            output.write("Codon\t")
            output.write("AA\t")
            output.write("Somatic VAF\t")
            output.write("Filter\t")
            output.write("COSMIC IDs\t")
            output.write("Num COSMIC Samples\t")
            output.write("COSMIC AA\t")
            output.write("Clinvar Significance\t")
            output.write("Clinvar HGVS\t")
            output.write("Clinvar Disease\t")
            output.write("Depth\t")
            output.write("Impact\t")
            output.write("Severity\t")
            output.write("Maximum Population AF\t")
            output.write("Chrom\t")
            output.write("Start\t")
            output.write("End\t")
            output.write("rsIDs\t")
            output.write("\n")
            for row in reader:
                sys.stdout.write("Processing sample {}\n".format(row[0]))

                vcfreader = cyvcf2.VCFReader("{}.vcfanno.vcf.gz".format(row[0]))
                desc = vcfreader["ANN"]["Description"]
                annotation_keys = [x.strip("\"'")
                                   for x in re.split("\s*\|\s*",
                                                     desc.split(":", 1)[1].strip('" '))]
                for v in VCF("{}.vcfanno.vcf.gz".format(row[0])):
                    if len(v.ALT) >= 1:
                        max_aaf = 1.0
                        if v.INFO.get('max_aaf_all'):
                            max_aaf = v.INFO.get('max_aaf_all')
                        if max_aaf < 0.005:
                            effects = utils.get_effects(v, annotation_keys)
                            top_impact = utils.get_top_impact(effects)
                            severity = top_impact.effect_severity
                            cosmic_data = utils.get_cosmic_info(v)
                            clinvar_data = get_clinvar_info(v)

                            output.write("{}\t".format(row[0]))
                            output.write("{}\t".format(top_impact.gene))
                            output.write("{}\t".format(v.REF))
                            output.write("{}\t".format(v.ALT))
                            output.write("{}\t".format(top_impact.codon_change))
                            output.write("{}\t".format(top_impact.aa_change))
                            output.write("{}\t".format(v.format('VF'))) # Somatic VAF
                            output.write("{}\t".format(v.FILTER))
                            output.write("{}\t".format(",".join(vcf_parsing.parse_cosmic_ids(variant))
                                                             or None))
                            output.write("{}\t".format(cosmic_data['num_samples']))
                            output.write("{}\t".format(cosmic_data['aa']))
                            output.write("{}\t".format(clinvar_data['significance']))
                            output.write("{}\t".format(clinvar_data['hgvs']))
                            output.write("{}\t".format(clinvar_data['disease']))
                            output.write("{}\t".format(v.INFO.get('DP')))
                            output.write("{}\t".format(top_impact.top_consequence))
                            output.write("{}\t".format(top_impact.effect_severity))
                            output.write("{}\t".format(max_aaf))
                            output.write("{}\t".format(v.CHROM))
                            output.write("{}\t".format(v.start))
                            output.write("{}\t".format(v.end))
                            output.write("\n")

    sys.stdout.write("Finished processing samples\n")

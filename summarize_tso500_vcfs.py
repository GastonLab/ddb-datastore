#!/usr/bin/env python

# This script is intended to act as a backup reporting system to convert VCF
# files directly to an excel report format, bypassing the VariantStore
# Cassandra database.

import re
import sys
import csv
import geneimpacts
import cyvcf2
import argparse

from cyvcf2 import VCF

def get_effects(variant, annotation_keys):
    effects = []
    effects += [geneimpacts.SnpEff(e, annotation_keys) for e in variant.INFO.get("ANN").split(",")]

    return effects


def get_top_impact(effects):
    top_impact = geneimpacts.Effect.top_severity(effects)

    if isinstance(top_impact, list):
        top_impact = top_impact[0]

    return top_impact

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

def get_cosmic_info(variant):
    cosmic_data = dict()

    cosmic_data['ids'] = variant.INFO.get('cosmic_ids') or 'None'
    cosmic_data['num_samples'] = variant.INFO.get('cosmic_numsamples') or 'None'
    cosmic_data['cds'] = variant.INFO.get('cosmic_cds') or 'None'
    cosmic_data['aa'] = variant.INFO.get('cosmic_aa') or 'None'
    cosmic_data['gene'] = variant.INFO.get('cosmic_gene') or 'None'

    return cosmic_data

def parse_cosmic_ids(variant_data):
    if variant_data.INFO.get('cosmic_ids') is not None:
        return variant_data.INFO.get('cosmic_ids').split(',')
    else:
        return []


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
                        if v.FILTER is None:
                            print("Variant passing filter found\n")
                            max_aaf = 1.0
                            if v.INFO.get('max_aaf_all'):
                                max_aaf = v.INFO.get('max_aaf_all')
                            if max_aaf < 0.005:
                                effects = get_effects(v, annotation_keys)
                                top_impact = get_top_impact(effects)
                                severity = top_impact.effect_severity
                                cosmic_data = get_cosmic_info(v)
                                clinvar_data = get_clinvar_info(v)

                                freq = v.format('VF')

                                output.write("{}\t".format(row[0]))
                                output.write("{}\t".format(top_impact.gene))
                                output.write("{}\t".format(v.REF))
                                output.write("{}\t".format(v.ALT))
                                output.write("{}\t".format(top_impact.codon_change))
                                output.write("{}\t".format(top_impact.aa_change))
                                output.write("{}\t".format(freq))
                                output.write("{}\t".format(v.FILTER))
                                output.write("{}\t".format(",".join(parse_cosmic_ids(v))
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

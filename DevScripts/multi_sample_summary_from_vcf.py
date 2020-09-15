#!/usr/bin/env python

# This script is intended to act as a backup reporting system to convert VCF
# files directly to an excel report format, bypassing the VariantStore
# Cassandra database and summarizing across samples from a run or runs

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


def process_sample_coverage(sample, samples, config):
    sample_coverage = defaultdict(lambda: dict())
    for library in samples[sample]:
        print sample
        print library
        report_panel_path = (
            "/mnt/shared-data/ddb-configs/disease_panels/{}/{}".format(
                samples[sample][library]['panel'],
                samples[sample][library]['report']))
        target_amplicons = utils.get_target_amplicons(report_panel_path)
        covf = "{}.sambamba_coverage.bed".format(samples[sample][library]['library_name'])
        with open(covf, 'rb') as coverage:
            reader = csv.reader(coverage, delimiter='\t')
            reader.next()

            for row in reader:
                if row[3] in target_amplicons:
                    sample_coverage[row[3]] = {
                        "num_reads": int(row[4]),
                        "mean_coverage": float(row[5])
                    }
    return sample_coverage


def process_sample_variants(sample, samples, config, thresholds):
    caller_records = defaultdict(lambda: dict())
    tier1_clinvar_terms = ("pathogenic", "likely-pathogenic", "drug-response")
    filtered_variant_data = defaultdict(list)

    for library in samples[sample]:
        report_panel_path = (
            "/mnt/shared-data/ddb-configs/disease_panels/{}/{}".format(
                samples[sample][library]['panel'],
                samples[sample][library]['report']))
        target_amplicons = utils.get_target_amplicons(report_panel_path)
        sys.stdout.write("Parsing Caller VCF Files\n")
        vcf_parsing.parse_vcf("{}.mutect.normalized.vcf.gz".format(samples[sample][library]['library_name']),
                              "mutect", caller_records)
        vcf_parsing.parse_vcf("{}.vardict.normalized.vcf.gz".format(samples[sample][library]['library_name']),
                              "vardict", caller_records)
        vcf_parsing.parse_vcf("{}.freebayes.normalized.vcf.gz".format(samples[sample][library]['library_name']),
                              "freebayes", caller_records)
        vcf_parsing.parse_vcf("{}.scalpel.normalized.vcf.gz".format(samples[sample][library]['library_name']),
                              "scalpel", caller_records)
        vcf_parsing.parse_vcf("{}.platypus.normalized.vcf.gz".format(samples[sample][library]['library_name']),
                              "platypus", caller_records)
        vcf_parsing.parse_vcf("{}.pindel.normalized.vcf.gz".format(samples[sample][library]['library_name']),
                              "pindel", caller_records)

        annotated_vcf = "{}.vcfanno.snpEff.GRCh37.75.vcf".format(samples[sample][library]['library_name'])

        sys.stdout.write("Parsing VCFAnno VCF\n")
        vcf = VCF(annotated_vcf)

        sys.stdout.write("Parsing VCFAnno VCF with CyVCF2\n")
        reader = cyvcf2.VCFReader(annotated_vcf)
        desc = reader["ANN"]["Description"]
        annotation_keys = [x.strip("\"'")
                           for x in re.split("\s*\|\s*",
                                             desc.split(":", 1)[1].strip('" '))]

        # Filter out variants with minor allele frequencies above the threshold but
        # retain any that are above the threshold but in COSMIC or in ClinVar and
        # not listed as benign.
        sys.stdout.write("Processing individual variants\n")
        for variant in vcf:
            if variant.INFO.get('max_aaf_all') < thresholds['max_maf']:
                callers = variant.INFO.get('CALLERS').split(',')
                effects = utils.get_effects(variant, annotation_keys)
                top_impact = utils.get_top_impact(effects)
                severity = top_impact.effect_severity
                amplicon_data = utils.get_amplicon_data(variant)
                amplicons = amplicon_data['amplicon'].split(',')

                assignable = 0
                for amplicon in amplicons:
                    if amplicon in target_amplicons:
                        assignable += 1
                        break
                if assignable:
                    key = (unicode("chr{}".format(variant.CHROM)),
                           int(variant.start), int(variant.end),
                           unicode(variant.REF), unicode(variant.ALT[0]))

                    caller_var_dicts = defaultdict(dict)
                    clinvar_data = utils.get_clinvar_info(variant, samples,
                                                          sample)
                    max_som_aaf = -1.00
                    max_depth = -1
                    min_depth = 100000000

                    for caller in callers:
                        caller_var_dicts[
                            caller] = parse_functions[caller](caller_records[caller][key])
                        if float(caller_var_dicts[caller]['AAF']) > max_som_aaf:
                            max_som_aaf = float(caller_var_dicts[caller]['AAF'])
                        if int(caller_var_dicts[caller]['DP']) < min_depth:
                            min_depth = int(caller_var_dicts[caller]['DP'])
                        if int(caller_var_dicts[caller]['DP']) > max_depth:
                            max_depth = int(caller_var_dicts[caller]['DP'])

                    if min_depth == 100000000:
                        min_depth = -1

                    # Putting in to Tier1 based on COSMIC
                    if vcf_parsing.var_is_in_cosmic(variant):
                        filtered_variant_data[
                                'tier1_variants'].append(variant)
                        continue
                    # Putting in to Tier1 based on ClinVar
                    if any(
                        i in tier1_clinvar_terms for i in clinvar_data[
                            'significance']):
                        filtered_variant_data[
                                'tier1_variants'].append(variant)
                        continue

                    if severity == 'MED' or severity == 'HIGH':
                        filtered_variant_data[
                                'tier3_variants'].append(variant)
                        continue
                    else:
                        filtered_variant_data[
                                'tier4_variants'].append(variant)
                        continue

    return filtered_variant_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file',
                        help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration',
                        help="Configuration file for various settings")
    parser.add_argument('-d', '--min_depth',
                        help='Minimum depth threshold for variant reporting',
                        default=200.0)
    parser.add_argument('-g', '--good_depth',
                        help='Floor for good depth of coverage',
                        default=500.0)
    parser.add_argument('-t', '--min_somatic_var_freq',
                        help='Minimum reportable somatic variant frequency',
                        default=0.01)
    parser.add_argument('-p', '--max_pop_freq',
                        help='Maximum allowed population allele frequency',
                        default=0.005)
    args = parser.parse_args()
    args.logLevel = "INFO"

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    libraries = configuration.configure_samples(args.samples_file, config)

    samples = configuration.merge_library_configs_samples(libraries)

    parse_functions = {'mutect': vcf_parsing.parse_mutect_vcf_record,
                       'freebayes': vcf_parsing.parse_freebayes_vcf_record,
                       'vardict': vcf_parsing.parse_vardict_vcf_record,
                       'scalpel': vcf_parsing.parse_scalpel_vcf_record,
                       'platypus': vcf_parsing.parse_platypus_vcf_record,
                       'pindel': vcf_parsing.parse_pindel_vcf_record}

    thresholds = {'min_saf': args.min_somatic_var_freq,
                  'max_maf': args.max_pop_freq,
                  'depth': args.min_depth}

    variant_summaries = defaultdict(lambda: defaultdict)
    coverage_summaries = defaultdict(lambda: defaultdict)

    for sample in samples:
        coverage_summaries[sample] = process_sample_coverage(sample, samples,
                                                             config)
        variant_summaries[sample] = process_sample_variants(sample, samples,
                                                            config, thresholds)

    sys.stdout.write("Finished processing samples\n")

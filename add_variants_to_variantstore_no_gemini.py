#!/usr/bin/env python

import re
import sys
import argparse

import utils
import cyvcf2
from cyvcf2 import VCF
from datetime import datetime
from collections import defaultdict
from cassandra.cqlengine import connection

from variantstore import Variant
from ddb import configuration
from ddb import vcf_parsing


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)")
    parser.add_argument('-v', '--variant_callers', help="Comma-delimited list of variant callers used")
    args = parser.parse_args()

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    connection.setup(['127.0.0.1'], "variantstore")

    parse_functions = {'mutect': vcf_parsing.parse_mutect_vcf_record,
                       'freebayes': vcf_parsing.parse_freebayes_vcf_record,
                       'vardict': vcf_parsing.parse_vardict_vcf_record,
                       'scalpel': vcf_parsing.parse_scalpel_vcf_record,
                       'platypus': vcf_parsing.parse_platypus_vcf_record,
                       'pindel': vcf_parsing.parse_pindel_vcf_record}

    thresholds = {'min_saf': 0.01,
                  'max_maf': 0.01,
                  'regions': config['actionable_regions']}

    for sample in samples:
        caller_vcf_records = defaultdict(lambda: dict())

        sys.stdout.write("Parsing Caller VCF Files")
        vcf_parsing.parse_vcf("{}.mutect.normalized.vcf".format(sample), "mutect", caller_vcf_records)
        vcf_parsing.parse_vcf("{}.vardict.normalized.vcf".format(sample), "vardict", caller_vcf_records)
        vcf_parsing.parse_vcf("{}.freebayes.normalized.vcf".format(sample), "freebayes", caller_vcf_records)
        vcf_parsing.parse_vcf("{}.scalpel.normalized.vcf".format(sample), "scalpel", caller_vcf_records)
        vcf_parsing.parse_vcf("{}.platypus.normalized.vcf".format(sample), "platypus", caller_vcf_records)
        vcf_parsing.parse_vcf("{}.pindel.normalized.vcf".format(sample), "pindel", caller_vcf_records)

        annotated_vcf = "{}.vcfanno.snpEff.GRCh37.75.vcf".format(sample)

        sys.stdout.write("Parsing VCFAnno VCF\n")
        vcf = VCF(annotated_vcf)

        sys.stdout.write("Parsing VCFAnno VCF with CyVCF2\n")
        reader = cyvcf2.VCFReader(annotated_vcf)
        desc = reader["ANN"]["Description"]
        annotation_keys = [x.strip("\"'") for x in re.split("\s*\|\s*", desc.split(":", 1)[1].strip('" '))]

        report_variants = list()

        # Filter out variants with minor allele frequencies above the threshold but
        # retain any that are above the threshold but in COSMIC or in ClinVar and not listed as benign.
        sys.stdout.write("Processing individual variants\n")
        for variant in vcf:
            effects = utils.get_effects(variant, annotation_keys)
            top_impact = utils.get_top_impact(effects)

            genes_list = utils.get_genes(effects)

            cassandra_variant = Variant(reference_genome=config['genome_version'],
                                        chr=variant.CHROM,
                                        pos=variant.start,
                                        end=variant.end,
                                        ref=variant.REF,
                                        alt=variant.ALT[0],
                                        sample=samples[sample]['sample_name'],
                                        extraction=samples[sample]['extraction'],
                                        library_name=sample,
                                        panel_name=samples[sample]['panel'],
                                        target_pool=samples[sample]['target_pool'],
                                        rs_id=variant.ID,
                                        date_annotated=datetime.now(),
                                        subtype=variant.INFO.get('sub_type'),
                                        type=variant.INFO.get('type'),
                                        gene=top_impact.gene,
                                        transcript=top_impact.transcript,
                                        exon=top_impact.exon,
                                        codon_change=top_impact.codon_change,
                                        biotype=top_impact.biotype,
                                        aa_change=top_impact.aa_change,
                                        severity=top_impact.effect_severity,
                                        impact=top_impact.top_consequence,
                                        impact_so=top_impact.so,
                                        max_maf_all=variant.INFO.get('max_aaf_all') or -1,
                                        max_maf_no_fin=variant.INFO.get('max_aaf_no_fin') or -1,
                                        # genes=utils.get_genes(effects),
                                        transcripts_data=utils.get_transcript_effects(effects),
                                        clinvar_data=utils.get_clinvar_info(variant),
                                        cosmic_data=utils.get_cosmic_info(variant),
                                        in_clinvar=vcf_parsing.var_is_in_clinvar(variant),
                                        in_cosmic=vcf_parsing.var_is_in_cosmic(variant),
                                        is_pathogenic=vcf_parsing.var_is_pathogenic(variant),
                                        is_lof=vcf_parsing.var_is_lof(variant),
                                        is_coding=vcf_parsing.var_is_coding(variant),
                                        is_splicing=vcf_parsing.var_is_splicing(variant),
                                        rs_ids=vcf_parsing.parse_rs_ids(variant),
                                        cosmic_ids=vcf_parsing.parse_cosmic_ids(variant)
                                        )

            if variant.INFO.get('CALLERS') is not None:
                cassandra_variant['callers'] = variant.INFO.get('CALLERS').split(',')

            population_freqs = {'esp_ea': variant.INFO.get('aaf_esp_ea') or -1,
                                'esp_aa': variant.INFO.get('aaf_esp_aa') or -1,
                                'esp_all': variant.INFO.get('aaf_esp_all') or -1,
                                '1kg_amr': variant.INFO.get('aaf_1kg_amr') or -1,
                                '1kg_eas': variant.INFO.get('aaf_1kg_eas') or -1,
                                '1kg_sas': variant.INFO.get('aaf_1kg_sas') or -1,
                                '1kg_afr': variant.INFO.get('aaf_1kg_afr') or -1,
                                '1kg_eur': variant.INFO.get('aaf_1kg_eur') or -1,
                                '1kg_all': variant.INFO.get('aaf_1kg_all') or -1,
                                'exac_all': variant.INFO.get('aaf_exac_all') or -1,
                                'adj_exac_all': variant.INFO.get('aaf_adj_exac_all') or -1,
                                'adj_exac_afr': variant.INFO.get('aaf_adj_exac_afr') or -1,
                                'adj_exac_amr': variant.INFO.get('aaf_adj_exac_amr') or -1,
                                'adj_exac_eas': variant.INFO.get('aaf_adj_exac_eas') or -1,
                                'adj_exac_fin': variant.INFO.get('aaf_adj_exac_fin') or -1,
                                'adj_exac_nfe': variant.INFO.get('aaf_adj_exac_nfe') or -1,
                                'adj_exac_oth': variant.INFO.get('aaf_adj_exac_oth') or -1,
                                'adj_exac_sas': variant.INFO.get('aaf_adj_exac_sas') or -1}

            cassandra_variant['population_freqs'] = population_freqs

            key = (unicode("chr{}".format(variant.CHROM)), int(variant.start), int(variant.end), unicode(variant.REF),
                   unicode(variant.ALT[0]))

            # sys.stdout.write("Processing data from individual callers\n")
            max_som_aaf = -1.00
            for caller in cassandra_variant['callers']:
                cassandra_variant[caller] = parse_functions[caller](caller_vcf_records[caller][key])
                if float(cassandra_variant[caller]['AAF']) > max_som_aaf:
                    max_som_aaf = float(cassandra_variant[caller]['AAF'])

            cassandra_variant['max_som_aaf'] = max_som_aaf

            report_variants.append(cassandra_variant)
            # sys.stdout.write("Saving data to Cassandra\n")
            cassandra_variant.save()

        if args.report:
            utils.write_sample_variant_report(args.report, sample, report_variants, args.variant_callers, thresholds)

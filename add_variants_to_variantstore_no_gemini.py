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
    args = parser.parse_args()

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    connection.setup(['127.0.0.1'], "variantstore")

    for sample in samples:
        caller_vcf_records = defaultdict(lambda: dict())

        vcf_parsing.parse_vcf(samples[sample]['mutect'], "mutect", caller_vcf_records)
        vcf_parsing.parse_vcf(samples[sample]['vardict'], "vardict", caller_vcf_records)
        vcf_parsing.parse_vcf(samples[sample]['freebayes'], "freebayes", caller_vcf_records)
        vcf_parsing.parse_vcf(samples[sample]['scalpel'], "scalpel", caller_vcf_records)

        vcf = VCF(samples[sample]['annotated_vcf'])

        reader = cyvcf2.VCFReader(samples[sample]['annotated_vcf'])
        desc = reader["ANN"]["Description"]
        annotation_keys = [x.strip("\"'") for x in re.split("\s*\|\s*", desc.split(":", 1)[1].strip('" '))]

        # Filter out variants with minor allele frequencies above the threshold but
        # retain any that are above the threshold but in COSMIC or in ClinVar and not listed as benign.
        for variant in vcf:
            effects = utils.get_effects(variant, annotation_keys)
            top_impact = utils.get_top_impact(effects)

            genes_list = utils.get_genes(effects)

            cassandra_variant = Variant(chr=variant.CHROM,
                                        start=variant.start,
                                        end=variant.end,
                                        ref=variant.REF,
                                        alt=variant.ALT[0],
                                        sample=samples[sample]['sample_name'],
                                        extraction=samples[sample]['extraction'],
                                        library_name=sample,
                                        panel_name=samples[sample]['panel'],
                                        target_pool=samples[sample]['target_pool'],
                                        rs_id=variant.ID,
                                        reference_genome=config['genome_version'],
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
                                        max_aaf_all=variant.INFO.get('max_aaf_all'),
                                        max_aaf_no_fin=variant.INFO.get('max_aaf_no_fin'),
                                        genes=utils.get_genes(effects),
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

            if 'mutect' in variant.INFO.get('CALLERS'):
                cassandra_variant['mutect'] = vcf_parsing.parse_mutect_vcf_record(caller_vcf_records['mutect'][key])

            if 'vardict' in variant.INFO.get('CALLERS'):
                cassandra_variant['vardict'] = vcf_parsing.parse_vardict_vcf_record(caller_vcf_records['vardict'][key])

            if 'freebayes' in variant.INFO.get('CALLERS'):
                cassandra_variant['freebayes'] = vcf_parsing.parse_vardict_vcf_record(caller_vcf_records['freebayes'][key])

            if 'scalpel' in variant.INFO.get('CALLERS'):
                cassandra_variant['scalpel'] = vcf_parsing.parse_vardict_vcf_record(caller_vcf_records['scalpel'][key])

            cassandra_variant.save()

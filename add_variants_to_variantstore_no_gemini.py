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
from ddb import gemini_interface
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
        annotation_parts = [x.strip("\"'") for x in re.split("\s*\|\s*", desc.split(":", 1)[1].strip('" '))]

        # Filter out variants with minor allele frequencies above the threshold but
        # retain any that are above the threshold but in COSMIC or in ClinVar and not listed as benign.
        for variant in vcf:
            impacts = utils.get_impacts(variant, annotation_parts)
            top_impact = utils.get_top_impact(impacts)

            cassandra_variant = Variant(chr=variant.INFO.get('chrom'),
                                        start=variant.INFO.get('start'),
                                        end=variant.INFO.get('end'),
                                        ref=variant.INFO.get('ref'),
                                        alt=variant.INFO.get('alt'),
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
                                        gene=top_impact['gene'],
                                        transcript=top_impact['transcript'],
                                        exon=top_impact['exon'],
                                        codon_change=top_impact['codon_change'],
                                        biotype=top_impact['biotype'],
                                        aa_change=top_impact['aa_change'],
                                        impact=top_impact['impact'],
                                        impact_so=top_impact['impact_so'],
                                        max_aaf_all=variant.INFO.get('max_aaf_all'),
                                        max_aaf_no_fin=variant.INFO.get('max_aaf_no_fin')
                                        )

            cassandra_variant['in_clinvar'] = gemini_interface.var_is_in_clinvar(variant.INFO.get)
            cassandra_variant['in_cosmic'] = gemini_interface.var_is_in_cosmic(variant.INFO.get)
            cassandra_variant['is_pathogenic'] = gemini_interface.var_is_pathogenic(variant.INFO.get)
            cassandra_variant['is_lof'] = gemini_interface.var_is_lof(variant.INFO.get)
            cassandra_variant['is_coding'] = gemini_interface.var_is_coding(variant.INFO.get)
            cassandra_variant['is_splicing'] = gemini_interface.var_is_splicing(variant.INFO.get)
            cassandra_variant['rs_ids'] = gemini_interface.parse_rs_ids(variant.INFO.get)
            cassandra_variant['cosmic_ids'] = gemini_interface.parse_cosmic_ids(variant.INFO.get)

            if variant.INFO.get('CALLERS') is not None:
                cassandra_variant['callers'] = variant.INFO.get('CALLERS').split(',')

            population_freqs = {'esp_ea': variant.INFO.get('aaf_esp_ea'),
                                'esp_aa': variant.INFO.get('aaf_esp_aa'),
                                'esp_all': variant.INFO.get('aaf_esp_all'),
                                '1kg_amr': variant.INFO.get('aaf_1kg_amr'),
                                '1kg_eas': variant.INFO.get('aaf_1kg_eas'),
                                '1kg_sas': variant.INFO.get('aaf_1kg_sas'),
                                '1kg_afr': variant.INFO.get('aaf_1kg_afr'),
                                '1kg_eur': variant.INFO.get('aaf_1kg_eur'),
                                '1kg_all': variant.INFO.get('aaf_1kg_all'),
                                'exac_all': variant.INFO.get('aaf_exac_all'),
                                'adj_exac_all': variant.INFO.get('aaf_adj_exac_all'),
                                'adj_exac_afr': variant.INFO.get('aaf_adj_exac_afr'),
                                'adj_exac_amr': variant.INFO.get('aaf_adj_exac_amr'),
                                'adj_exac_eas': variant.INFO.get('aaf_adj_exac_eas'),
                                'adj_exac_fin': variant.INFO.get('aaf_adj_exac_fin'),
                                'adj_exac_nfe': variant.INFO.get('aaf_adj_exac_nfe'),
                                'adj_exac_oth': variant.INFO.get('aaf_adj_exac_oth'),
                                'adj_exac_sas': variant.INFO.get('aaf_adj_exac_sas')}

            cassandra_variant['population_freqs'] = population_freqs

            key = (variant.INFO.get('chrom'), variant.INFO.get('start'), variant.INFO.get('end'), 
                   variant.INFO.get('ref'), variant.INFO.get('alt'))

            if 'mutect' in variant.INFO.get('CALLERS'):
                cassandra_variant['mutect'] = vcf_parsing.parse_mutect_vcf_record(caller_vcf_records['mutect'][key])

            if 'vardict' in variant.INFO.get('CALLERS'):
                cassandra_variant['vardict'] = vcf_parsing.parse_vardict_vcf_record(caller_vcf_records['vardict'][key])

            if 'freebayes' in variant.INFO.get('CALLERS'):
                cassandra_variant['freebayes'] = vcf_parsing.parse_vardict_vcf_record(caller_vcf_records['freebayes'][key])

            if 'scalpel' in variant.INFO.get('CALLERS'):
                cassandra_variant['scalpel'] = vcf_parsing.parse_vardict_vcf_record(caller_vcf_records['scalpel'][key])

            cassandra_variant.save()

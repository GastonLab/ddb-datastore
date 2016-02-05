#!/usr/bin/env python

import sys
import argparse
import configuration

from datetime import datetime
from gemini import GeminiQuery
from collections import defaultdict
from cassandra.cqlengine import connection

from variantstore import Variant
from cyvcf2 import VCF


def _var_is_rare(variant_data):
    """Determine if the MAF of the variant is < 1% in all populations"""

    if variant_data['in_esp'] != 0 or variant_data['in_1kg'] != 0 or variant_data['in_exac'] != 0:
        if variant_data['max_aaf_all'] > 0.01:
            return False
        else:
            return True
    else:
        return True


def _var_is_in_cosmic(variant_data):
    """Determine if the variant has a COSMIC identifier"""

    if variant_data['cosmic_ids'] is not None:
        return True
    else:
        return False


def _var_is_in_clinvar(variant_data):
    """Determine if there is ClinVar data for the variant"""

    if variant_data['clinvar_sig'] is not None:
        return True
    else:
        return False


def _var_is_pathogenic(variant_data):
    """Determine if a variant in clinvar is potentially pathogenic"""

    if variant_data['clinvar_sig'] is not None:
        if "pathogenic" in variant_data['clinvar_sig']:
            return True
        else:
            return False
    else:
        return False


def _var_is_protein_effecting(variant_data):
    if variant_data['impact_severity'] != "LOW":
        return True
    else:
        return False


def _var_in_gene(variant_data, genes):
    if variant_data['gene'] in genes:
        return True
    else:
        return False


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
        query = "SELECT chrom, start, end, ref, alt, vcf_id, rs_ids, cosmic_ids, filter, qual, qual_depth, depth, " \
                "type, sub_type, " \
                "gene, transcript, exon, codon_change, aa_change, biotype, impact, impact_so, impact_severity, " \
                "aa_length, is_lof, is_conserved, pfam_domain, in_omim, clinvar_sig, clinvar_disease_name, " \
                "is_exonic, is_coding, is_splicing, " \
                "clinvar_origin, clinvar_causal_allele, clinvar_dbsource, clinvar_dbsource_id, " \
                "clinvar_on_diag_assay, rmsk, in_segdup, strand_bias, rms_map_qual, in_hom_run, num_mapq_zero, " \
                "num_reads_w_dels, grc, gms_illumina, in_cse, num_alleles, allele_count, haplotype_score, " \
                "is_somatic, somatic_score, aaf_esp_ea, aaf_esp_aa, aaf_esp_all, aaf_1kg_amr, " \
                "aaf_1kg_eas, aaf_1kg_sas, aaf_1kg_afr, aaf_1kg_eur, aaf_1kg_all, aaf_exac_all, aaf_adj_exac_all, " \
                "aaf_adj_exac_afr, aaf_adj_exac_amr, aaf_adj_exac_eas, aaf_adj_exac_fin, aaf_adj_exac_nfe, " \
                "aaf_adj_exac_oth, aaf_adj_exac_sas, max_aaf_all, in_esp, in_1kg, in_exac, info," \
                "(gts).(*), (gt_depths).(*), (gt_ref_depths).(*), (gt_alt_depths).(*) FROM variants"

        sys.stdout.write("Running GEMINI query\n")
        gq = GeminiQuery(samples[sample]['db'])
        gq.run(query)
        header = gq.header

        caller_vcf_records = defaultdict(lambda: dict())

        sys.stdout.write("Reading {}\n".format(samples[sample]['mutect']))
        mutect_vcf = VCF(samples[sample]['mutect'])
        for record in mutect_vcf:
            if len(record.ALT) > 1:
                sys.stderr.write("ERROR: More than one alternative allele detected in file "
                                 "{}\n Record: {}\n".format(samples[sample]['mutect'], record))
                sys.exit()
            key = (unicode("chr{}".format(record.CHROM)), int(record.start), int(record.end), unicode(record.REF),
                   unicode(record.ALT[0]))
            caller_vcf_records['mutect'][key] = record

        sys.stdout.write("Reading {}\n".format(samples[sample]['vardict']))
        vardict_vcf = VCF(samples[sample]['vardict'])
        for record in vardict_vcf:
            if len(record.ALT) > 1:
                sys.stderr.write("ERROR: More than one alternative allele detected in file "
                                 "{}\n Record: {}\n".format(samples[sample]['varddict'], record))
                sys.exit()
            key = (unicode("chr{}".format(record.CHROM)), int(record.start), int(record.end), unicode(record.REF),
                   unicode(record.ALT[0]))
            caller_vcf_records['vardict'][key] = record

        sys.stdout.write("Reading {}\n".format(samples[sample]['freebayes']))
        freebayes_vcf = VCF(samples[sample]['freebayes'])
        for record in freebayes_vcf:
            if len(record.ALT) > 1:
                sys.stderr.write("ERROR: More than one alternative allele detected in file "
                                 "{}\n Record: {}\n".format(samples[sample]['freebayes'], record))
                sys.exit()
            key = (unicode("chr{}".format(record.CHROM)), int(record.start), int(record.end), unicode(record.REF),
                   unicode(record.ALT[0]))
            caller_vcf_records['freebayes'][key] = record

        sys.stdout.write("Reading {}\n".format(samples[sample]['scalpel']))
        scalpel_vcf = VCF(samples[sample]['scalpel'])
        for record in scalpel_vcf:
            if len(record.ALT) > 1:
                sys.stderr.write("ERROR: More than one alternative allele detected in file "
                                 "{}\n Record: {}\n".format(samples[sample]['scalpel'], record))
                sys.exit()
            key = (unicode("chr{}".format(record.CHROM)), int(record.start), int(record.end), unicode(record.REF),
                   unicode(record.ALT[0]))
            caller_vcf_records['scalpel'][key] = record

        # Filter out variants with minor allele frequencies above the threshold but
        # retain any that are above the threshold but in COSMIC or in ClinVar and not listed as benign.
        for variant_data in gq:
            cassandra_variant = Variant(chr=variant_data['chrom'], start=variant_data['start'], end=variant_data['end'],
                                        ref=variant_data['ref'], alt=variant_data['alt'],
                                        sample=samples[sample]['sample_name'], extraction=samples[sample]['extraction'],
                                        library_name=sample, panel_name=samples[sample]['panel'],
                                        target_pool=samples[sample]['target_pool'],
                                        reference_genome=config['genome_version'], date_annotated=datetime.now(),
                                        type=variant_data['type'],
                                        subtype=variant_data['sub_type'], rs_id=variant_data['vcf_id'],
                                        gene=variant_data['gene'], max_aaf=variant_data['max_aaf_all'],
                                        transcript=variant_data['transcript'], exon=variant_data['exon'],
                                        codon_change=variant_data['codon_change'],
                                        aa_change=variant_data['aa_change'], biotype=variant_data['biotype'],
                                        impact=variant_data['impact'], impact_so=variant_data['impact_so'])

            cassandra_variant['in_clinvar'] = _var_is_in_clinvar(variant_data)
            cassandra_variant['in_cosmic'] = _var_is_in_cosmic(variant_data)
            cassandra_variant['is_pathogenic'] = _var_is_pathogenic(variant_data)

            if variant_data['is_lof']:
                cassandra_variant['is_lof'] = True
            else:
                cassandra_variant['is_lof'] = False

            if variant_data['is_coding']:
                cassandra_variant['is_coding'] = True
            else:
                cassandra_variant['is_coding'] = False

            if variant_data['is_splicing']:
                cassandra_variant['is_splicing'] = True
            else:
                cassandra_variant['is_splicing'] = False

            if variant_data['rs_ids'] is not None:
                cassandra_variant['rs_ids'] = variant_data['rs_ids'].split(',')

            if variant_data['cosmic_ids'] is not None:
                cassandra_variant['cosmic_ids'] = variant_data['cosmic_ids'].split(',')

            if variant_data['info']['CALLERS'] is not None:
                cassandra_variant['callers'] = variant_data['info']['CALLERS'].split(',')

            population_freqs = {'esp_ea': variant_data['aaf_esp_ea'],
                                'esp_aa': variant_data['aaf_esp_aa'],
                                'esp_all': variant_data['aaf_esp_all'],
                                '1kg_amr': variant_data['aaf_1kg_amr'],
                                '1kg_eas': variant_data['aaf_1kg_eas'],
                                '1kg_sas': variant_data['aaf_1kg_sas'],
                                '1kg_afr': variant_data['aaf_1kg_afr'],
                                '1kg_eur': variant_data['aaf_1kg_eur'],
                                '1kg_all': variant_data['aaf_1kg_all'],
                                'exac_all': variant_data['aaf_exac_all'],
                                'adj_exac_all': variant_data['aaf_adj_exac_all'],
                                'adj_exac_afr': variant_data['aaf_adj_exac_afr'],
                                'adj_exac_amr': variant_data['aaf_adj_exac_amr'],
                                'adj_exac_eas': variant_data['aaf_adj_exac_eas'],
                                'adj_exac_fin': variant_data['aaf_adj_exac_fin'],
                                'adj_exac_nfe': variant_data['aaf_adj_exac_nfe'],
                                'adj_exac_oth': variant_data['aaf_adj_exac_oth'],
                                'adj_exac_sas': variant_data['aaf_adj_exac_sas']}

            cassandra_variant['population_freqs'] = population_freqs

            key = (variant_data['chrom'], variant_data['start'], variant_data['end'], variant_data['ref'],
                   variant_data['alt'])

            if 'mutect' in variant_data['info']['CALLERS']:
                record = caller_vcf_records['mutect'][key]
                # filters = record.FILTER.split(',')

                mutect_info = {'GTF_DP': str(record.gt_depths[0]),
                               'GTF_AD': str(record.gt_alt_depths[0])}

                cassandra_variant['mutect'] = mutect_info

            if 'vardict' in variant_data['info']['CALLERS']:
                record = caller_vcf_records['vardict'][key]
                # filters = record.FILTER.split(',')

                vardict_info = {'DP': str(record.INFO.get('DP')),
                                'VD': str(record.INFO.get('VD')),
                                'AF': str(record.INFO.get('AF')),
                                'BIAS': str(record.INFO.get('BIAS')),
                                'REFBIAS': str(record.INFO.get('REFBIAS')),
                                'VARBIAS': str(record.INFO.get('VARBIAS')),
                                'QUAL': str(record.INFO.get('QUAL')),
                                'QSTD': str(record.INFO.get('QSTD')),
                                'SBF': str(record.INFO.get('SBF')),
                                'ODDRATIO': str(record.INFO.get('ODDRATIO')),
                                'MQ': str(record.INFO.get('MQ')),
                                'SN': str(record.INFO.get('SN')),
                                'HIAF': str(record.INFO.get('HIAF')),
                                'ADJAF': str(record.INFO.get('ADJAF')),
                                'MSI': str(record.INFO.get('MSI')),
                                'MSILEN': str(record.INFO.get('MSILEN')),
                                'SHIFT3': str(record.INFO.get('SHIFT3')),
                                'NM': str(record.INFO.get('NM')),
                                'GDAMP': str(record.INFO.get('GDAMP')),
                                'LSEQ': str(record.INFO.get('LSEQ')),
                                'RSEQ': str(record.INFO.get('RSEQ')),
                                'TLAMP': str(record.INFO.get('TLAMP')),
                                'NCAMP': str(record.INFO.get('NCAMP')),
                                'AMPFLAG': str(record.INFO.get('AMPFLAG')),
                                'HICNT': str(record.INFO.get('HICNT')),
                                'HICOV': str(record.INFO.get('HICOV')),
                                'GTF_DP': str(record.gt_depths[0]),
                                'GTF_AD': str(record.gt_alt_depths[0])}

                cassandra_variant['vardict'] = vardict_info

            if 'freebayes' in variant_data['info']['CALLERS']:
                record = caller_vcf_records['freebayes'][key]
                # filters = record.FILTER.split(',')

                freebayes_info = {'DP': str(record.INFO.get('DP')),
                                  'AF': str(record.INFO.get('AF')),
                                  'AC': str(record.INFO.get('AC')),
                                  'RO': str(record.INFO.get('RO')),
                                  'AO': str(record.INFO.get('AO')),
                                  'PRO': str(record.INFO.get('PRO')),
                                  'PAO': str(record.INFO.get('PAO')),
                                  'QR': str(record.INFO.get('QR')),
                                  'QA': str(record.INFO.get('QA')),
                                  'PQR': str(record.INFO.get('PQR')),
                                  'PQA': str(record.INFO.get('PQA')),
                                  'SRF': str(record.INFO.get('SRF')),
                                  'SRR': str(record.INFO.get('SRR')),
                                  'SAF': str(record.INFO.get('SAF')),
                                  'SAR': str(record.INFO.get('SAR')),
                                  'SRP': str(record.INFO.get('SRP')),
                                  'SAP': str(record.INFO.get('SAP')),
                                  'AB': str(record.INFO.get('AB')),
                                  'ABP': str(record.INFO.get('ABP')),
                                  'RUN': str(record.INFO.get('RUN')),
                                  'RPP': str(record.INFO.get('RPP')),
                                  'RPPR': str(record.INFO.get('RPPR')),
                                  'RPL': str(record.INFO.get('RPL')),
                                  'RPR': str(record.INFO.get('RPR')),
                                  'EPP': str(record.INFO.get('EPP')),
                                  'EPPR': str(record.INFO.get('EPPR')),
                                  'DRPA': str(record.INFO.get('DRPA')),
                                  'ODDS': str(record.INFO.get('ODDS')),
                                  'GTI': str(record.INFO.get('GTI')),
                                  'TYPE': str(record.INFO.get('TYPE')),
                                  'CIGAR': str(record.INFO.get('CIGAR')),
                                  'NUMALT': str(record.INFO.get('NUMALT')),
                                  'MEANALT': str(record.INFO.get('MEANALT')),
                                  'LEN': str(record.INFO.get('LEN')),
                                  'MQM': str(record.INFO.get('MQM')),
                                  'MQMR': str(record.INFO.get('MQMR')),
                                  'PAIRED': str(record.INFO.get('PAIRED')),
                                  'PAIREDR': str(record.INFO.get('PAIREDR')),
                                  'GTF_DP': str(record.gt_depths[0])}

                cassandra_variant['freebayes'] = freebayes_info

            if 'scalpel' in variant_data['info']['CALLERS']:
                record = caller_vcf_records['scalpel'][key]
                # filters = record.FILTER.split(',')

                scalpel_info = {'AVGCOV': str(record.INFO.get('AVGCOV')),
                                'MINCOV': str(record.INFO.get('MINCOV')),
                                'ALTCOV': str(record.INFO.get('ALTCOV')),
                                'COVRATIO': str(record.INFO.get('COVRATIO')),
                                'ZYG': str(record.INFO.get('ZYG')),
                                'CHI2': str(record.INFO.get('CHI2')),
                                'FISHERPHREDSCORE': str(record.INFO.get('FISHERPHREDSCORE')),
                                'INH': str(record.INFO.get('INH')),
                                'BESTSTATE': str(record.INFO.get('BESTSTATE')),
                                'COVSTATE': str(record.INFO.get('COVSTATE')),
                                'SOMATIC': str(record.INFO.get('SOMATIC')),
                                'DENOVO': str(record.INFO.get('DENOVO')),
                                'GTF_DP': str(record.gt_depths[0]),
                                'GTF_AD': str(record.gt_alt_depths[0])}

                cassandra_variant['scalpel'] = scalpel_info

            # for key in cassandra_variant.keys():
            #     sys.stdout.write("{}: {} ({})\n".format(key, cassandra_variant[key], type(cassandra_variant[key])))
            cassandra_variant.save()

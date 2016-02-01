#!/usr/bin/env python

import sys
import cyvcf
import argparse
import ConfigParser
from gemini import GeminiQuery
from collections import defaultdict
from cassandra.cqlengine import connection

from variantstore import Variant
from cyvcf2 import VCF


def configure_runtime(infile):
    """Parse the configuration settings from a file
    :param infile: input filename
    :type infile: string.
    :returns:  dict -- A configuration dictionary.
    """

    configuration = defaultdict()
    config = ConfigParser.SafeConfigParser()
    config.read(infile)

    try:
        config.options('settings')
    except ConfigParser.NoSectionError:
        sys.stderr.write("No section: settings in file\n")
        sys.exit()

    try:
        config.options('resources')
    except ConfigParser.NoSectionError:
        sys.stderr.write("No section: resources in file\n")
        sys.exit()

    # Set all options specified in file
    for option in config.options('settings'):
            configuration[option] = config.get('settings', option)

    for resource in config.options('resources'):
            configuration[resource] = config.get('resources', resource)

    # Configure individual tools
    for section in config.sections():
        if section != 'settings' and section != 'resources':
            tool = section
            options = config.options(tool)
            tool_dict = dict()

            # Set all specified options
            for option in options:
                tool_dict[option] = config.get(tool, option)

            configuration[tool] = tool_dict

    return configuration


def configure_samples(infile, configuration):
    """Parse the sample-level configuration settings from a file
    :param infile: input filename
    :type infile: string.
    :param configuration: project/run level configuration
    :type configuration: dictionary
    :returns:  dict -- A configuration dictionary.
    """

    samples = dict()

    config = ConfigParser.SafeConfigParser()
    config.read(infile)

    for sample in config.sections():
        sample_dict = dict()
        for option in config.options(sample):
            sample_dict[option] = config.get(sample, option)

        if 'regions' not in sample_dict.keys():
            sample_dict['regions'] = configuration['regions']
        if 'snv_regions' not in sample_dict.keys():
            if 'snv_regions' in configuration:
                sample_dict['snv_regions'] = configuration['snv_regions']
        if 'indel_regions' not in sample_dict.keys():
            if 'indel_regions' in configuration:
                sample_dict['indel_regions'] = configuration['indel_regions']

        samples[sample] = sample_dict

    return samples


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    args = parser.parse_args()

    sys.stdout.write("Parsing configuration data\n")
    config = configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configure_samples(args.samples_file, config)

    for sample in samples:
        query = "SELECT chrom, start, end, ref, alt, vcf_id, rs_ids, cosmic_ids, filter, qual, qual_depth, depth, " \
                "type, sub_type, " \
                "gene, transcript, exon, codon_change, aa_change, biotype, impact, impact_so, impact_severity, " \
                "aa_length, is_lof, is_conserved, pfam_domain, in_omim, clinvar_sig, clinvar_disease_name, " \
                "clinvar_origin, clinvar_causal_allele, clinvar_dbsource, clinvar_dbsource_id, " \
                "clinvar_on_diag_assay, rmsk, in_segdup, strand_bias, rms_map_qual, in_hom_run, num_mapq_zero, " \
                "num_reads_w_dels, grc, gms_illumina, in_cse, num_alleles, allele_count, haplotype_score, " \
                "is_somatic, somatic_score, aaf_esp_ea, aaf_esp_aa, aaf_esp_aa, aaf_esp_all, aaf_1kg_amr, " \
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
                                        ref=variant_data['ref'], alt=variant_data['alt'], sample=sample,
                                        callers=variant_data['info']['CALLERS'], type=variant_data['type'],
                                        subtype=variant_data['sub_type'], rs_id=variant_data['vcf_id'],
                                        rs_ids=variant_data['rs_ids'].split(','),
                                        cosmic_ids=variant_data['cosmic_ids'].split(','), gene=variant_data['gene'],
                                        transcript=variant_data['transcript'], exon=variant_data['exon'],
                                        codon_change=variant_data['codon_change'],
                                        aa_change=variant_data['aa_change'], biotype=variant_data['biotype'],
                                        impact=variant_data['impact'], impact_so=variant_data['impact_so'],
                                        max_aaf=variant_data['max_aaf'])

            key = (variant_data['chrom'], variant_data['start'], variant_data['end'], variant_data['ref'],
                   variant_data['alt'])

            if 'mutect' in variant_data['info']['CALLERS']:
                caller_info['mutect'] = caller_vcf_records['mutect'][key]
                sys.stdout.write("MuTect: ")
                print caller_vcf_records['mutect'][key]

            if 'vardict' in variant_data['info']['CALLERS']:
                caller_info['vardict'] = caller_vcf_records['vardict'][key]
                sys.stdout.write("VarDict: ")
                print caller_vcf_records['vardict'][key]

            if 'freebayes' in variant_data['info']['CALLERS']:
                caller_info['freebayes'] = caller_vcf_records['freebayes'][key]
                sys.stdout.write("FreeBayes: ")
                print caller_vcf_records['freebayes'][key]

            if 'scalpel' in variant_data['info']['CALLERS']:
                caller_info['scalpel'] = caller_vcf_records['scalpel'][key]
                sys.stdout.write("Scalpel: ")
                print caller_vcf_records['scalpel'][key]

            sys.stdout.write("******************************************\n")
            cassandra_variant.save()

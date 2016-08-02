#!/usr/bin/env python

import re
import sys
import getpass
import argparse
import utils
import cyvcf2

from cyvcf2 import VCF
from datetime import datetime
from collections import defaultdict
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider

from variantstore import Variant
from variantstore import SampleVariant
from variantstore import TargetVariant
from ddb import configuration
from ddb import vcf_parsing
from ddb_ngsflow import pipeline

from toil.job import Job


def process_sample(job, addresses, keyspace, authenticator, parse_functions, variant_callers,
                   sample, samples, config):
    connection.setup(addresses, keyspace, auth_provider=authenticator)

    caller_records = defaultdict(lambda: dict())

    sys.stdout.write("Parsing Caller VCF Files\n")
    vcf_parsing.parse_vcf("{}.mutect.normalized.vcf".format(sample), "mutect", caller_records)
    vcf_parsing.parse_vcf("{}.vardict.normalized.vcf".format(sample), "vardict", caller_records)
    vcf_parsing.parse_vcf("{}.freebayes.normalized.vcf".format(sample), "freebayes", caller_records)
    vcf_parsing.parse_vcf("{}.scalpel.normalized.vcf".format(sample), "scalpel", caller_records)
    vcf_parsing.parse_vcf("{}.platypus.normalized.vcf".format(sample), "platypus", caller_records)
    vcf_parsing.parse_vcf("{}.pindel.normalized.vcf".format(sample), "pindel", caller_records)

    annotated_vcf = "{}.vcfanno.snpEff.GRCh37.75.vcf".format(sample)

    sys.stdout.write("Parsing VCFAnno VCF\n")
    vcf = VCF(annotated_vcf)

    sys.stdout.write("Parsing VCFAnno VCF with CyVCF2\n")
    reader = cyvcf2.VCFReader(annotated_vcf)
    desc = reader["ANN"]["Description"]
    annotation_keys = [x.strip("\"'") for x in re.split("\s*\|\s*", desc.split(":", 1)[1].strip('" '))]

    # Filter out variants with minor allele frequencies above the threshold but
    # retain any that are above the threshold but in COSMIC or in ClinVar and not listed as benign.
    sys.stdout.write("Processing individual variants\n")
    for variant in vcf:
        # Parsing VCF and creating data structures for Cassandra model
        callers = variant.INFO.get('CALLERS').split(',')
        effects = utils.get_effects(variant, annotation_keys)
        top_impact = utils.get_top_impact(effects)
        population_freqs = get_population_freqs(variant)
        amplicon_data = get_amplicon_data(variant)

        key = (unicode("chr{}".format(variant.CHROM)), int(variant.start), int(variant.end), unicode(variant.REF),
               unicode(variant.ALT[0]))

        caller_variant_data_dicts = defaultdict(dict)
        max_som_aaf = -1.00
        max_depth = -1
        min_depth = 100000000

        for caller in callers:
            caller_variant_data_dicts[caller] = parse_functions[caller](caller_records[caller][key])
            if float(caller_variant_data_dicts[caller]['AAF']) > max_som_aaf:
                max_som_aaf = float(caller_variant_data_dicts[caller]['AAF'])
            if int(caller_variant_data_dicts[caller]['DP']) < min_depth:
                min_depth = int(caller_variant_data_dicts[caller]['DP'])
            if int(caller_variant_data_dicts[caller]['DP']) > max_depth:
                max_depth = int(caller_variant_data_dicts[caller]['DP'])

        if min_depth == 100000000:
            min_depth = -1

        # Create Cassandra Objects
        # Create the general variant ordered table
        cassandra_variant = Variant.create(
                reference_genome=config['genome_version'],
                chr=variant.CHROM,
                pos=variant.start,
                end=variant.end,
                ref=variant.REF,
                alt=variant.ALT[0],
                sample=samples[sample]['sample_name'],
                extraction=samples[sample]['extraction'],
                library_name=samples[sample]['library_name'],
                run_id=samples[sample]['run_id'],
                panel_name=samples[sample]['panel'],
                target_pool=samples[sample]['target_pool'],
                sequencer=samples[sample]['sequencer'],
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
                cosmic_ids=vcf_parsing.parse_cosmic_ids(variant),
                callers=callers,
                population_freqs=population_freqs,
                amplicon_data=amplicon_data,
                max_som_aaf=max_som_aaf,
                min_depth=min_depth,
                max_depth=max_depth,
                mutect=caller_variant_data_dicts['mutect'] or dict(),
                freebayes=caller_variant_data_dicts['freebayes'] or dict(),
                scalpel=caller_variant_data_dicts['scalpel'] or dict(),
                platypus=caller_variant_data_dicts['platypus'] or dict(),
                pindel=caller_variant_data_dicts['pindel'] or dict(),
                vardict=caller_variant_data_dicts['vardict'] or dict(),
                manta=caller_variant_data_dicts['manta'] or dict()
                )

        # Create Cassandra Object
        sample_variant = SampleVariant.create(
                sample=samples[sample]['sample_name'],
                run_id=samples[sample]['run_id'],
                library_name=samples[sample]['library_name'],
                reference_genome=config['genome_version'],
                chr=variant.CHROM,
                pos=variant.start,
                end=variant.end,
                ref=variant.REF,
                alt=variant.ALT[0],
                extraction=samples[sample]['extraction'],
                panel_name=samples[sample]['panel'],
                target_pool=samples[sample]['target_pool'],
                sequencer=samples[sample]['sequencer'],
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
                cosmic_ids=vcf_parsing.parse_cosmic_ids(variant),
                callers=callers,
                population_freqs=population_freqs,
                amplicon_data=amplicon_data,
                max_som_aaf=max_som_aaf,
                min_depth=min_depth,
                max_depth=max_depth,
                mutect=caller_variant_data_dicts['mutect'] or dict(),
                freebayes=caller_variant_data_dicts['freebayes'] or dict(),
                scalpel=caller_variant_data_dicts['scalpel'] or dict(),
                platypus=caller_variant_data_dicts['platypus'] or dict(),
                pindel=caller_variant_data_dicts['pindel'] or dict(),
                vardict=caller_variant_data_dicts['vardict'] or dict(),
                manta=caller_variant_data_dicts['manta'] or dict()
                )

        target_variant = TargetVariant.create(
            target=amplicon_data['amplicon'],
            sample=samples[sample]['sample_name'],
            run_id=samples[sample]['run_id'],
            library_name=samples[sample]['library_name'],
            reference_genome=config['genome_version'],
            chr=variant.CHROM,
            pos=variant.start,
            end=variant.end,
            ref=variant.REF,
            alt=variant.ALT[0],
            extraction=samples[sample]['extraction'],
            panel_name=samples[sample]['panel'],
            target_pool=samples[sample]['target_pool'],
            sequencer=samples[sample]['sequencer'],
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
            cosmic_ids=vcf_parsing.parse_cosmic_ids(variant),
            callers=callers,
            population_freqs=population_freqs,
            amplicon_data=amplicon_data,
            max_som_aaf=max_som_aaf,
            min_depth=min_depth,
            max_depth=max_depth,
            mutect=caller_variant_data_dicts['mutect'] or dict(),
            freebayes=caller_variant_data_dicts['freebayes'] or dict(),
            scalpel=caller_variant_data_dicts['scalpel'] or dict(),
            platypus=caller_variant_data_dicts['platypus'] or dict(),
            pindel=caller_variant_data_dicts['pindel'] or dict(),
            vardict=caller_variant_data_dicts['vardict'] or dict(),
            manta=caller_variant_data_dicts['manta'] or dict()
        )

    job.fileStore.logToMaster("Data saved to Cassandra for sample {}\n".format(sample))


def get_population_freqs(variant):
    freqs = {'esp_ea': variant.INFO.get('aaf_esp_ea') or -1,
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

    return freqs


def get_amplicon_data(variant):
    data = {'amplicon': variant.INFO.get('amplicon_target') or "None",
            'intersect': variant.INFO.get('amplicon_intersect') or "None"}

    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-v', '--variant_callers', help="Comma-delimited list of variant callers used")
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    args.logLevel = "INFO"

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
    else:
        auth_provider = None

    parse_functions = {'mutect': vcf_parsing.parse_mutect_vcf_record,
                       'freebayes': vcf_parsing.parse_freebayes_vcf_record,
                       'vardict': vcf_parsing.parse_vardict_vcf_record,
                       'scalpel': vcf_parsing.parse_scalpel_vcf_record,
                       'platypus': vcf_parsing.parse_platypus_vcf_record,
                       'pindel': vcf_parsing.parse_pindel_vcf_record}

    root_job = Job.wrapJobFn(pipeline.spawn_batch_jobs, cores=1)

    for sample in samples:
        sample_job = Job.wrapJobFn(process_sample, [args.address], "variantstore", auth_provider, parse_functions,
                                   args.variant_callers, sample, samples, config,
                                   cores=1)
        root_job.addChild(sample_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)

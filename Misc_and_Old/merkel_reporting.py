#!/usr/bin/env python

import sys
import utils
import argparse
import getpass

from ddb import configuration
from collections import defaultdict

from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)", default='any_evidence_report')
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)

    parser.add_argument('-d', '--min_depth', help='Minimum depth threshold for variant reporting', default=100.0)
    parser.add_argument('-t', '--min_somatic_var_freq', help='Minimum reportable somatic variant frequency',
                        default=0.10)
    parser.add_argument('-p', '--max_pop_freq', help='Maximum allowed population allele frequency', default=0.005)

    args = parser.parse_args()
    args.logLevel = "INFO"

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    libraries = configuration.configure_samples(args.samples_file, config)
    samples = configuration.merge_library_configs_samples(libraries)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
        connection.setup([args.address], "variantstore", auth_provider=auth_provider)
    else:
        connection.setup([args.address], "variantstore")

    thresholds = {'min_saf': float(args.min_somatic_var_freq),
                  'max_maf': float(args.max_pop_freq),
                  'depth': float(args.min_depth)}

    callers = ("mutect", "platypus", "vardict", "scalpel", "freebayes", "pindel")
    project_variant_data = defaultdict(lambda: defaultdict(int))
    variant_count_data = defaultdict(lambda: defaultdict(int))
    gene_count_data = defaultdict(lambda: defaultdict(int))
    variants_list = list()

    sys.stdout.write("Processing samples\n")
    for sample in samples:
        sys.stdout.write("Processing variants for sample {}\n".format(sample))
        target_amplicon_coverage = defaultdict(lambda: defaultdict(float))

        report_names = {'log': "{}.{}.log".format(sample, args.report),
                        'coverage': "{}_coverage_{}.txt".format(sample, args.report),
                        'tier1_pass': "{}_tier1_pass_variants_{}.txt".format(sample, args.report),
                        'vus_pass': "{}_vus_pass_variants_{}.txt".format(sample, args.report),
                        'tier4_pass': "{}_tier4_pass_variants_{}.txt".format(sample, args.report),
                        'all_ordered': "{}_all_ordered_variants_{}.txt".format(sample, args.report),
                        'categories': "{}_variants_by_category{}.txt".format(sample, args.report)
                        }

        with open(report_names['log'], 'w') as logfile:
            logfile.write("Reporting Log for sample {}\n".format(sample))
            logfile.write("---------------------------------------------\n")

        with open(report_names['coverage'], "w") as coverage_report:
            coverage_report.write("Sample:\t{}\n".format(sample))
            coverage_report.write("---------------------------------------------\n")

        utils.setup_report_header(report_names['tier1_pass'], callers)
        utils.setup_report_header(report_names['vus_pass'], callers)
        utils.setup_report_header(report_names['tier4_pass'], callers)
        utils.setup_report_header(report_names['all_ordered'], callers)

        for library in samples[sample]:
            report_panel_path = "/mnt/shared-data/ddb-configs/disease_panels/{}/{}" \
                                "".format(samples[sample][library]['panel'], samples[sample][library]['report'])
            target_amplicons = utils.get_target_amplicons(report_panel_path)

            with open(report_names['log'], 'a') as logfile:
                sys.stdout.write("Processing amplicons for library {} from file {}\n".format(library,
                                                                                             report_panel_path))
                logfile.write("Processing amplicons for library {} from file {}\n".format(library, report_panel_path))

            ordered_variants, num_var = utils.get_variants(config, samples, sample, library, thresholds, report_names)

            sys.stdout.write("Processing amplicon coverage data\n")
            reportable_amplicons, target_amplicon_coverage = utils.get_coverage_data(target_amplicons, samples, sample,
                                                                                     library, target_amplicon_coverage)

            sys.stdout.write("Filtering and classifying variants\n")
            filtered_var_data = utils.classify_and_filter_variants_proj(samples, sample, library, report_names,
                                                                        target_amplicons, callers, ordered_variants,
                                                                        config, thresholds, project_variant_data,
                                                                        variant_count_data, gene_count_data,
                                                                        variants_list, args.address, auth_provider)
            variants_list = filtered_var_data[-4]
            project_variant_data = filtered_var_data[-3]
            variant_count_data = filtered_var_data[-2]
            gene_count_data = filtered_var_data[-1]

            sys.stdout.write("Writing variant reports\n")
            utils.write_reports(report_names, samples, sample, library, filtered_var_data, ordered_variants,
                                target_amplicon_coverage, reportable_amplicons, num_var, thresholds, callers)

    sys.stdout.write("Writing project/run level data\n")
    with open("Summary_Data.txt", 'w') as summary:
        summary.write("Variant\tNum Tier1 Pass\tNum Tier1 Fail\tNum VUS Pass\tNum VUS Fail\tNum Tier4 Pass\t"
                      "Num Tier4 Fail\n")
        for variant_id in project_variant_data:
            summary.write("{}\t{}\t{}\t{}\n".format(variant_id,
                                                    project_variant_data[variant_id]['tier1_pass'],
                                                    project_variant_data[variant_id]['vus_pass'],
                                                    project_variant_data[variant_id]['tier4_pass']))

    sys.stdout.write("Writing project/run level category data\n")
    with open("Category_Data.txt", 'w') as summary:
        summary.write("Variant\tNum Pos\tNum Neg\tDiff\n")
        for variant_id in project_variant_data:
            diff = abs(project_variant_data[variant_id]['positive'] - project_variant_data[variant_id]['negative'])
            summary.write("{}\t{}\t{}\t{}\n".format(variant_id, project_variant_data[variant_id]['positive'],
                                                    project_variant_data[variant_id]['negative'], diff))

    sys.stdout.write("Writing Sample-level variant count data\n")
    with open("Sample_Variant_Counts.txt", 'w') as summary:
        summary.write("sample\tgroup\tviral_status\tnum_pass\tnum_ct\tnum_high_pathogenic\tnum_med\tnum_low\n")
        for sample in variant_count_data:
            for library in samples[sample]:
                summary.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, samples[sample][library]['category'],
                                                                samples[sample][library]['viral'],
                                                                variant_count_data[sample]['pass_count'],
                                                                variant_count_data[sample]['CT_count'],
                                                                variant_count_data[sample]['high_impact_pathogenic'],
                                                                variant_count_data[sample]['med'],
                                                                variant_count_data[sample]['low']))

    sys.stdout.write("Writing Sample and gene-level variant count data\n")
    with open("Sample_Gene_Variant_Counts.txt", 'w') as summary:
        summary.write("Sample\tGroup\tViral Status\tRB1\tTP53\tNOTCH2\tMLL3\tBRCA1\tPIK3CA\tAPC\tNOTCH1\tMLL2\n")
        for sample in gene_count_data:
            for library in samples[sample]:
                summary.write("{}\t{}\t{}\t{}\t{}\n".format(sample, samples[sample][library]['category'],
                                                            samples[sample][library]['viral'],
                                                            gene_count_data[sample]['RB1'],
                                                            gene_count_data[sample]['TP53'],
                                                            gene_count_data[sample]['NOTCH2'],
                                                            gene_count_data[sample]['MLL3'],
                                                            gene_count_data[sample]['BRCA1'],
                                                            gene_count_data[sample]['PIK3CA'],
                                                            gene_count_data[sample]['APC'],
                                                            gene_count_data[sample]['NOTCH1'],
                                                            gene_count_data[sample]['MLL2']))

    sys.stdout.write("Writing variant data\n")
    with open("Variant_Data.txt", 'w') as summary:
        summary.write("Sample\tCategory\tViral Status\tGene\tAmplicon\tRef\tAlt\tCodon\tAA\tImpact\tSeverity\t"
                      "Fraction\tMax AF\tCallers\tCOSMIC IDs\tCOSMIC Num Samples\tCOSMIC AA\tClinVar Sig\tHGVS\t"
                      "COSMIC Disease\tMean Cov\tNum Reads\tMax MAF\tMin Depth\tMax Depth\tChr\tStart\tEnd\trsIDs\n"
                      "\n".format())
        for variant in variants_list:
            summary.write("{sample}\t{cat}\t{status}\t{gene}\t{amp}\t{ref}\t{alt}\t{codon}\t{aa}\t{impact}\t"
                          "{severity}\t{frac}\t{max_som_aaf}\t{callers}\t{cosmic}\t{cosmic_nsamples}\t{cosmic_aa}\t"
                          "{csig}\t{hgvs}\t{cdis}\t{cov}\t{reads}\t{max_maf_all}\t{min_depth}\t{max_depth}\t{chr}\t"
                          "{start}\t{end}\t{rsids}\n"
                          "".format(sample=variant.sample,
                                    cat=samples[variant.sample][variant.library_name]['category'],
                                    status=samples[variant.sample][variant.library_name]['viral'],
                                    chr=variant.chr,
                                    start=variant.pos,
                                    end=variant.end,
                                    gene=variant.gene,
                                    ref=variant.ref,
                                    alt=variant.alt,
                                    codon=variant.codon_change,
                                    aa=variant.aa_change,
                                    rsids=",".join(variant.rs_ids),
                                    cosmic=",".join(variant.cosmic_ids) or None,
                                    cosmic_nsamples=variant.cosmic_data['num_samples'],
                                    cosmic_aa=variant.cosmic_data['aa'],
                                    amp=variant.amplicon_data['amplicon'],
                                    csig=variant.clinvar_data['significance'],
                                    hgvs=variant.clinvar_data['hgvs'],
                                    cdis=variant.clinvar_data['disease'],
                                    cov=target_amplicon_coverage[variant.amplicon_data['amplicon']]['mean_coverage'],
                                    reads=target_amplicon_coverage[variant.amplicon_data['amplicon']]['num_reads'],
                                    impact=variant.impact,
                                    frac=variant.fraction,
                                    severity=variant.severity,
                                    max_maf_all=variant.max_maf_all,
                                    max_som_aaf=variant.max_som_aaf,
                                    min_depth=variant.min_depth,
                                    max_depth=variant.max_depth,
                                    callers=",".join(variant.callers) or None))



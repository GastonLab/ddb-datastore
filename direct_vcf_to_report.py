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
                        "num_reads": row[4],
                        "mean_coverage": row[5]
                    }
    return sample_coverage


def process_sample_variants(coverage, sample, samples, config, thresholds):
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
                        if max_som_aaf < thresholds['min_saf']:
                            filtered_variant_data[
                                'tier1_fail_variants'].append(variant)
                        elif max_depth < thresholds['depth']:
                            filtered_variant_data[
                                'tier1_fail_variants'].append(variant)
                        else:
                            filtered_variant_data[
                                'tier1_pass_variants'].append(variant)
                        continue
                    # Putting in to Tier1 based on ClinVar
                    if any(
                        i in tier1_clinvar_terms for i in clinvar_data[
                            'significance']):
                        if max_som_aaf < thresholds['min_saf']:
                            filtered_variant_data[
                                'tier1_fail_variants'].append(variant)
                        elif max_depth < thresholds['depth']:
                            filtered_variant_data[
                                'tier1_fail_variants'].append(variant)
                        else:
                            filtered_variant_data[
                                'tier1_pass_variants'].append(variant)
                        continue

                    if severity == 'MED' or severity == 'HIGH':
                        if max_som_aaf < thresholds['min_saf']:
                            filtered_variant_data[
                                'tier3_fail_variants'].append(variant)
                        elif max_depth < thresholds['depth']:
                            filtered_variant_data[
                                'tier3_fail_variants'].append(variant)
                        else:
                            filtered_variant_data[
                                'tier3_pass_variants'].append(variant)
                        continue
                    else:
                        if max_som_aaf < thresholds['min_saf']:
                            filtered_variant_data[
                                'tier4_fail_variants'].append(variant)
                        elif max_depth < thresholds['depth']:
                            filtered_variant_data[
                                'tier4_fail_variants'].append(variant)
                        else:
                            filtered_variant_data[
                                'tier4_pass_variants'].append(variant)
                        continue

    report_name = "{}.xlsx".format(sample)

    wb = xlwt.Workbook()

    error_style = xlwt.easyxf('pattern: pattern solid, fore_colour red;')
    warning_style = xlwt.easyxf(
        'pattern: pattern solid, fore_colour light_orange;')
    pass_style = xlwt.easyxf(
        'pattern: pattern solid, fore_colour light_green;')
    default_style = xlwt.easyxf('pattern: pattern solid, fore_colour white;')

    coverage_sheet = wb.add_sheet("Coverage")
    tier1_sheet = wb.add_sheet("Tier1 and 2 Pass")
    tier1_fail_sheet = wb.add_sheet("Tier1 and 2 Fail")
    tier3_sheet = wb.add_sheet("Tier3 Pass")
    tier3_fail_sheet = wb.add_sheet("Tier3 Fail")
    tier4_sheet = wb.add_sheet("Tier4 Pass")
    tier4_fail_sheet = wb.add_sheet("Tier4 Fail")

    tier_sheets = (tier1_sheet, tier1_fail_sheet, tier3_sheet,
                   tier3_fail_sheet, tier4_sheet, tier4_fail_sheet)
    tier_key = ("tier1_pass_variants", "tier1_fail_variants",
                "tier3_pass_variants", "tier3_fail_variants",
                "tier4_pass_variants", "tier4_fail_variants")

    libraries = list()
    report_templates = list()
    run_id = ""
    for library in samples[sample]:
        libraries.append(samples[sample][library]['library_name'])
        report_templates.append(samples[sample][library]['report'])
        run_id = samples[sample][library]['run_id']
    lib_string = " | ".join(libraries)
    reports_string = " | ".join(report_templates)

    coverage_sheet.write(0, 0, "Sample")
    coverage_sheet.write(0, 1, "{}".format(sample))

    coverage_sheet.write(1, 0, "Libraries")
    coverage_sheet.write(1, 1, "{}".format(lib_string))

    coverage_sheet.write(2, 0, "Run ID")
    coverage_sheet.write(2, 1, "{}".format(run_id))

    coverage_sheet.write(3, 0, "Reporting Templates")
    coverage_sheet.write(3, 1, "{}".format(reports_string))

    coverage_sheet.write(4, 0, "Minimum Reportable Somatic Allele Frequency")
    coverage_sheet.write(4, 1, "{}".format(thresholds['min_saf']))

    coverage_sheet.write(5, 0, "Minimum Amplicon Depth")
    coverage_sheet.write(5, 1, "{}".format(thresholds['depth']))

    coverage_sheet.write(6, 0, "Maximum Population Allele Frequency")
    coverage_sheet.write(6, 1, "{}".format(thresholds['max_maf']))

    coverage_sheet.write(7, 0, "Amplicon")
    coverage_sheet.write(7, 1, "Num Reads")
    coverage_sheet.write(7, 2, "Coverage")

    row_num = 8
    for amplicon in target_amplicons:
        print amplicon
        if coverage[amplicon]['mean_coverage'] < 200:
            style = error_style
        elif coverage[amplicon]['mean_coverage'] < 500:
            style = warning_style
        else:
            style = pass_style

        coverage_sheet.write(row_num, 0, "{}".format(amplicon),
                             style)
        coverage_sheet.write(row_num, 1, "{}".format(coverage[amplicon
                                                              ]['num_reads']),
                             style)
        coverage_sheet.write(row_num, 2, "{}".format(coverage[amplicon
                                                              ]['mean_coverage']),
                             style)

        row_num += 1

    ###########################################################################

    sheet_num = 0
    for sheet in tier_sheets:
        sheet.write(0, 0, "Gene")
        sheet.write(0, 1, "Amplicon")
        sheet.write(0, 2, "Ref")
        sheet.write(0, 3, "Alt")
        sheet.write(0, 4, "Codon")
        sheet.write(0, 5, "AA")
        sheet.write(0, 6, "Max Caller Somatic VAF")
        sheet.write(0, 7, "Callers")
        sheet.write(0, 8, "Caller Counts")
        sheet.write(0, 9, "COSMIC IDs")
        sheet.write(0, 10, "Num COSMIC Samples")
        sheet.write(0, 11, "COSMIC AA")
        sheet.write(0, 12, "Clinvar Significance")
        sheet.write(0, 13, "Clinvar HGVS")
        sheet.write(0, 14, "Clinvar Disease")
        sheet.write(0, 15, "Coverage")
        sheet.write(0, 16, "Num Reads")
        sheet.write(0, 17, "Impact")
        sheet.write(0, 18, "Severity")
        sheet.write(0, 19, "Maximum Population AF")
        sheet.write(0, 20, "Min Caller Depth")
        sheet.write(0, 21, "Max Caller Depth")
        sheet.write(0, 22, "Chrom")
        sheet.write(0, 23, "Start")
        sheet.write(0, 24, "End")
        sheet.write(0, 25, "rsIDs")

        col = 26
        if 'mutect' in callers:
            sheet.write(0, col, "MuTect_AF")
            col += 1

        if 'vardict' in callers:
            sheet.write(0, col, "VarDict_AF")
            col += 1

        if 'freebayes' in callers:
            sheet.write(0, col, "FreeBayes_AF")
            col += 1

        if 'scalpel' in callers:
            sheet.write(0, col, "Scalpel_AF")
            col += 1

        if 'platypus' in callers:
            sheet.write(0, col, "Platypus_AF")
            col += 1

        if 'pindel' in callers:
            sheet.write(0, col, "Pindel_AF")
            col += 1

        row = 1
        for variant in filtered_variant_data[tier_key[sheet_num]]:
            callers = variant.INFO.get('CALLERS').split(',')
            effects = utils.get_effects(variant, annotation_keys)
            top_impact = utils.get_top_impact(effects)
            severity = top_impact.effect_severity
            amplicon_data = utils.get_amplicon_data(variant)
            amplicons = amplicon_data['amplicon'].split(',')
            clinvar_data = utils.get_clinvar_info(variant, samples, sample)
            cosmic_data = utils.get_cosmic_info(variant)
            max_aaf_all = variant.INFO.get('max_aaf_all') or -1

            if "pathogenic" in clinvar_data['significance']:
                style = pass_style
            elif "drug-response" in clinvar_data['significance']:
                style = pass_style
            elif "likely-pathogenic" in clinvar_data['significance']:
                style = pass_style
            else:
                style = default_style

            coverage_values = list()
            reads_values = list()
            for amplicon in amplicons:
                coverage_values.append(
                    str(coverage[amplicon]['mean_coverage']))
                reads_values.append(
                    str(coverage[amplicon]['num_reads']))

            coverage_string = ",".join(coverage_values)
            reads_string = ",".join(reads_values)

            if len(variant.ref) < 200:
                ref = variant.REF
            else:
                ref = "Length > 200bp"

            if len(variant.alt) < 200:
                alt = variant.ALT[0]
            else:
                alt = "Length > 200bp"

            if len(top_impact.codon_change) < 200:
                codon_change = top_impact.codon_change
            else:
                codon_change = "Length > 200aa"

            if len(top_impact.aa_change) < 200:
                aa_change = top_impact.aa_change
            else:
                aa_change = "Length > 200aa"

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

            sheet.write(row, 0, "{}".format(top_impact.gene), style)
            sheet.write(row, 1, "{}".format(amplicon_data['amplicon']),
                        style)
            sheet.write(row, 2, "{}".format(ref), style)
            sheet.write(row, 3, "{}".format(alt), style)
            sheet.write(row, 4, "{}".format(codon_change), style)
            sheet.write(row, 5, "{}".format(aa_change), style)
            sheet.write(row, 6, "{}".format(max_som_aaf), style)
            sheet.write(row, 7, "{}".format(variant.INFO.get('CALLERS')
                                            or None), style)
            sheet.write(row, 8, "{}".format(variant.num_times_callers), style)
            sheet.write(row, 9, "{}".format(",".join(vcf_parsing.parse_cosmic_ids(variant))
                                             or None), style)
            sheet.write(row, 10,
                        "{}".format(cosmic_data['num_samples']), style)
            sheet.write(row, 11, "{}".format(cosmic_data['aa']), style)
            sheet.write(row, 12,
                        "{}".format(clinvar_data['significance']),
                        style)
            sheet.write(row, 13,
                        "{}".format(clinvar_data['hgvs']), style)
            sheet.write(row, 14,
                        "{}".format(clinvar_data['disease']), style)
            sheet.write(row, 15, "{}".format(coverage_string), style)
            sheet.write(row, 16, "{}".format(reads_string), style)
            sheet.write(row, 17, "{}".format(top_impact.top_consequence),
                        style)
            sheet.write(row, 18, "{}".format(top_impact.effect_severity),
                        style)
            sheet.write(row, 19, "{}".format(max_aaf_all), style)
            sheet.write(row, 20, "{}".format(min_depth), style)
            sheet.write(row, 21, "{}".format(max_depth), style)
            sheet.write(row, 22, "{}".format(variant.CHROM), style)
            sheet.write(row, 23, "{}".format(variant.start), style)
            sheet.write(row, 24, "{}".format(variant.end), style)
            sheet.write(row, 25, "{}".format(",".join(
                vcf_parsing.parse_rs_ids(variant))), style)

            col = 26
            if 'mutect' in callers:
                sheet.write(row, col, "{}".format(variant.mutect.get('AAF')
                                                  or None), style)
                col += 1

            if 'vardict' in callers:
                sheet.write(row, col, "{}".format(variant.vardict.get('AAF')
                                                  or None), style)
                col += 1

            if 'freebayes' in callers:
                sheet.write(row, col, "{}".format(variant.freebayes.get('AAF')
                                                  or None), style)
                col += 1

            if 'scalpel' in callers:
                sheet.write(row, col, "{}".format(variant.scalpel.get('AAF')
                                                  or None), style)
                col += 1

            if 'platypus' in callers:
                sheet.write(row, col, "{}".format(variant.platypus.get('AAF')
                                                  or None), style)
                col += 1

            if 'pindel' in callers:
                sheet.write(row, col, "{}".format(variant.pindel.get('AAF')
                                                  or None), style)
                col += 1

            row += 1
        sheet_num += 1
    wb.save(report_name)


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

    for sample in samples:
        sample_coverage = process_sample_coverage(sample, samples, config)
        process_sample_variants(sample_coverage, sample, samples, config,
                                thresholds)

    sys.stdout.write("Finished processing samples\n")

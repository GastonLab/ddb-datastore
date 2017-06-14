#!/usr/bin/env python

import sys
import xlwt
import utils
import getpass
import argparse

import numpy as np

from toil.job import Job
from ddb import configuration
from ddb_ngsflow import pipeline
from variantstore import Variant
from collections import defaultdict
from variantstore import SampleVariant
from coveragestore import SampleCoverage
from cassandra.cqlengine import connection
from cassandra.auth import PlainTextAuthProvider


def process_sample(job, config, sample, samples, addresses, authenticator, thresholds, callers):
    job.fileStore.logToMaster("Retrieving data for sample {}\n".format(sample))
    job.fileStore.logToMaster("Retrieving coverage data from database\n")
    connection.setup(addresses, "coveragestore", auth_provider=authenticator)

    report_data = dict()
    filtered_variant_data = defaultdict(list)
    off_target_amplicon_counts = defaultdict(int)
    target_amplicon_coverage = dict()
    ordered_amplicon_coverage = list()

    iterated = 0
    passing_variants = 0
    filtered_low_freq = 0
    filtered_low_depth = 0
    filtered_off_target = 0

    for library in samples[sample]:
        counted = list()
        category = samples[sample][library]['category']
        report_panel_path = "/mnt/shared-data/ddb-configs/disease_panels/{}/{}" \
                            "".format(samples[sample][library]['panel'], samples[sample][library]['report'])

        job.fileStore.logToMaster("{}: processing amplicons from file {}".format(library, report_panel_path))
        target_amplicons = utils.get_target_amplicons(report_panel_path)

        for amplicon in target_amplicons:
            coverage_data = SampleCoverage.objects.timeout(None).filter(
                SampleCoverage.sample == samples[sample][library]['sample_name'],
                SampleCoverage.amplicon == amplicon,
                SampleCoverage.run_id == samples[sample][library]['run_id'],
                SampleCoverage.library_name == samples[sample][library]['library_name'],
                SampleCoverage.program_name == "sambamba"
            )
            ordered_amplicons = coverage_data.order_by('amplicon', 'run_id').limit(coverage_data.count() + 1000)
            for result in ordered_amplicons:
                target_amplicon_coverage[amplicon] = result
                ordered_amplicon_coverage.append(result)

        job.fileStore.logToMaster("{}: retrieving variants".format(library))
        variants = SampleVariant.objects.timeout(None).filter(
            SampleVariant.reference_genome == config['genome_version'],
            SampleVariant.sample == samples[sample][library]['sample_name'],
            SampleVariant.run_id == samples[sample][library]['run_id'],
            SampleVariant.library_name == samples[sample][library]['library_name'],
            SampleVariant.max_maf_all <= thresholds['max_maf']
        ).allow_filtering()

        num_var = variants.count()
        ordered = variants.order_by('library_name', 'chr', 'pos', 'ref', 'alt').limit(variants.count() + 1000)
        job.fileStore.logToMaster("{}: retrieved {} variants from database\n".format(library, num_var))
        job.fileStore.logToMaster("{}: classifying and filtering variants\n".format(library))

        for variant in ordered:
            iterated += 1
            if len(variant.callers) < 2:
                continue
            if len(variant.ref) > 2 and len(variant.alt) > 2:
                continue

            elements = variant.amplicon_data['amplicon'].split('_')
            gene = elements[0]
            variant_id = "{}:{}-{}_{}_{}_{}_{}".format(variant.chr, variant.pos, variant.end, variant.ref, variant.alt,
                                                       variant.codon_change, variant.aa_change)

            if variant.amplicon_data['amplicon'] is 'None':
                filtered_off_target += 1
                off_target_amplicon_counts[variant.amplicon_data['amplicon']] += 1
            else:
                amplicons = variant.amplicon_data['amplicon'].split(',')
                assignable = 0
                for amplicon in amplicons:
                    if amplicon in target_amplicons:
                        assignable += 1
                        break
                if assignable:
                    if variant.max_som_aaf > thresholds['min_saf']:
                        if variant.min_depth > thresholds['depth']:
                            if variant_id not in counted:

                                match_variants = Variant.objects.timeout(None).filter(
                                    Variant.reference_genome == config['genome_version'],
                                    Variant.chr == variant.chr,
                                    Variant.pos == variant.pos,
                                    Variant.ref == variant.ref,
                                    Variant.alt == variant.alt
                                ).allow_filtering()

                                num_matches = match_variants.count()
                                ordered_var = match_variants.order_by('ref', 'alt', 'sample', 'library_name',
                                                                      'run_id').limit(num_matches + 1000)
                                vafs = list()
                                num_times_callers = defaultdict(int)
                                for var in ordered_var:
                                    vaf = var.max_som_aaf
                                    vafs.append(vaf)
                                    for caller in var.callers:
                                        num_times_callers[caller] += 1

                                variant.vaf_median = np.median(vafs)
                                variant.vaf_std_dev = np.std(vafs)
                                variant.num_times_called = num_matches


                                counted.append(variant_id)

                                caller_counts_elements = list()
                                for caller in num_times_callers:
                                    caller_counts_elements.append("{}: {}".format(caller, num_times_callers[caller]))
                                variant.num_times_callers = ",".join(caller_counts_elements)

                                # Putting in to Tier1 based on COSMIC
                                if variant.cosmic_ids:
                                    if variant.max_som_aaf < thresholds['min_saf']:
                                        filtered_variant_data['tier1_fail_variants'].append(variant)
                                        filtered_low_freq += 1
                                    elif variant.max_depth < thresholds['depth']:
                                        filtered_variant_data['tier1_fail_variants'].append(variant)
                                        filtered_low_depth += 1
                                    else:
                                        filtered_variant_data['tier1_pass_variants'].append(variant)
                                        passing_variants += 1
                                    continue

                                # Putting in to Tier1 based on ClinVar not being None or Benign
                                if variant.clinvar_data['pathogenic'] != 'None':
                                    if variant.clinvar_data['pathogenic'] != 'benign':
                                        if variant.clinvar_data['pathogenic'] != 'likely-benign':
                                            if variant.max_som_aaf < thresholds['min_saf']:
                                                filtered_variant_data['tier1_fail_variants'].append(variant)
                                                filtered_low_freq += 1
                                            elif variant.max_depth < thresholds['depth']:
                                                filtered_variant_data['tier1_fail_variants'].append(variant)
                                                filtered_low_depth += 1
                                            else:
                                                filtered_variant_data['tier1_pass_variants'].append(variant)
                                                passing_variants += 1
                                            continue

                                if variant.severity == 'MED' or variant.severity == 'HIGH':
                                    if variant.max_som_aaf < thresholds['min_saf']:
                                        filtered_variant_data['tier3_fail_variants'].append(variant)
                                        filtered_low_freq += 1
                                    elif variant.max_depth < thresholds['depth']:
                                        filtered_variant_data['tier3_fail_variants'].append(variant)
                                        filtered_low_depth += 1
                                    else:
                                        filtered_variant_data['tier3_pass_variants'].append(variant)
                                        passing_variants += 1
                                    continue
                                else:
                                    if variant.max_som_aaf < thresholds['min_saf']:
                                        filtered_variant_data['tier4_fail_variants'].append(variant)
                                        filtered_low_freq += 1
                                    elif variant.max_depth < thresholds['depth']:
                                        filtered_variant_data['tier4_fail_variants'].append(variant)
                                        filtered_low_depth += 1
                                    else:
                                        filtered_variant_data['tier4_pass_variants'].append(variant)
                                        passing_variants += 1
                                    continue
                else:
                    filtered_off_target += 1
                    off_target_amplicon_counts[variant.amplicon_data['amplicon']] += 1

        job.fileStore.logToMaster("{}: iterated through {} variants\n".format(library, iterated))
        job.fileStore.logToMaster("{}: filtered {} off-target variants\n".format(library, filtered_off_target))
        job.fileStore.logToMaster("{}: filtered {} low-freq variants\n".format(library, filtered_low_freq))
        job.fileStore.logToMaster("{}: filtered {} low-depth variants\n".format(library, filtered_low_depth))

        job.fileStore.logToMaster("{}: passing {} tier 1 and 2 variants"
                                  "\n".format(library, len(filtered_variant_data['tier1_pass_variants'])))
        job.fileStore.logToMaster("{}: passing {} tier3 variants"
                                  "\n".format(library, len(filtered_variant_data['tier3_pass_variants'])))
        job.fileStore.logToMaster("{}: passing {} tier 4 variants"
                                  "\n".format(library, len(filtered_variant_data['tier4_pass_variants'])))

    report_data['variants'] = filtered_variant_data
    report_data['coverage'] = target_amplicon_coverage

    report_name = "{}.xlsx".format(sample)

    wb = xlwt.Workbook()

    error_style = xlwt.easyxf('pattern: pattern solid, fore_colour red;')
    warning_style = xlwt.easyxf('pattern: pattern solid, fore_colour light_orange;')
    pass_style = xlwt.easyxf('pattern: pattern solid, fore_colour light_green;')

    coverage_sheet = wb.add_sheet("Coverage")
    tier1_sheet = wb.add_sheet("Tier1 and 2 Pass")
    tier3_sheet = wb.add_sheet("Tier3 Pass")
    tier4_sheet = wb.add_sheet("Tier4 Pass")
    tier1_fail_sheet = wb.add_sheet("Tier1 and 2 Fail")
    tier3_fail_sheet = wb.add_sheet("Tier3 Fail")
    tier4_fail_sheet = wb.add_sheet("Tier4 Fail")

    tier_sheets = (tier1_sheet, tier1_fail_sheet, tier3_sheet, tier3_fail_sheet, tier4_sheet, tier4_fail_sheet)
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

    coverage_sheet.write(7, 0, "Sample")
    coverage_sheet.write(7, 1, "Library")
    coverage_sheet.write(7, 2, "Amplicon")
    coverage_sheet.write(7, 3, "Num Reads")
    coverage_sheet.write(7, 4, "Coverage")

    row_num = 9
    for amplicon in ordered_amplicon_coverage:
        if amplicon.mean_coverage < 250:
            style = error_style
        elif amplicon.mean_coverage < 500:
            style = warning_style
        else:
            style = pass_style

        coverage_sheet.write(row_num, 0, "{}".format(amplicon.sample), style)
        coverage_sheet.write(row_num, 1, "{}".format(amplicon.library_name), style)
        coverage_sheet.write(row_num, 2, "{}".format(amplicon.amplicon), style)
        coverage_sheet.write(row_num, 3, "{}".format(amplicon.num_reads), style)
        coverage_sheet.write(row_num, 4, "{}".format(amplicon.mean_coverage), style)

        row_num += 1

    ####################################################################################################################

    sheet_num = 0
    for sheet in tier_sheets:
        sheet.write(0, 0, "Sample")
        sheet.write(0, 1, "Library")
        sheet.write(0, 2, "Gene")
        sheet.write(0, 3, "Amplicon")
        sheet.write(0, 4, "Ref")
        sheet.write(0, 5, "Alt")
        sheet.write(0, 6, "Codon")
        sheet.write(0, 7, "AA")
        sheet.write(0, 8, "Max Caller Somatic VAF")
        sheet.write(0, 9, "Num Times in Database")
        sheet.write(0, 10, "Median VAF")
        sheet.write(0, 11, "StdDev VAF")
        sheet.write(0, 12, "Callers")
        sheet.write(0, 13, "Caller Counts")
        sheet.write(0, 14, "COSMIC IDs")
        sheet.write(0, 15, "Num COSMIC Samples")
        sheet.write(0, 16, "COSMIC AA")
        sheet.write(0, 17, "Clinvar Significance")
        sheet.write(0, 18, "Clinvar HGVS")
        sheet.write(0, 19, "Clinvar Disease")
        sheet.write(0, 20, "Coverage")
        sheet.write(0, 21, "Num Reads")
        sheet.write(0, 22, "Impact")
        sheet.write(0, 23, "Severity")
        sheet.write(0, 24, "Maximum Population AF")
        sheet.write(0, 25, "Min Caller Depth")
        sheet.write(0, 26, "Max Caller Depth")
        sheet.write(0, 27, "Chrom")
        sheet.write(0, 28, "Start")
        sheet.write(0, 29, "End")
        sheet.write(0, 30, "rsIDs")

        col = 31
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
        for variant in report_data['variants'][tier_key[sheet_num]]:
            amplicons = variant.amplicon_data['amplicon'].split(',')

            coverage_values = list()
            reads_values = list()
            for amplicon in amplicons:
                coverage_values.append(str(report_data['coverage'][amplicon]['mean_coverage']))
                reads_values.append(str(report_data['coverage'][amplicon]['num_reads']))

            coverage_string = ",".join(coverage_values)
            reads_string = ",".join(reads_values)

            sheet.write(row, 0, "{}".format(variant.sample))
            sheet.write(row, 1, "{}".format(variant.library_name))
            sheet.write(row, 2, "{}".format(variant.gene))
            sheet.write(row, 3, "{}".format(variant.amplicon_data['amplicon']))
            sheet.write(row, 4, "{}".format(variant.ref))
            sheet.write(row, 5, "{}".format(variant.alt))
            sheet.write(row, 6, "{}".format(variant.codon_change))
            sheet.write(row, 7, "{}".format(variant.aa_change))
            sheet.write(row, 8, "{}".format(variant.max_som_aaf))
            sheet.write(row, 9, "{}".format(variant.num_times_called))
            sheet.write(row, 10, "{}".format(variant.vaf_median))
            sheet.write(row, 11, "{}".format(variant.vaf_std_dev))
            sheet.write(row, 12, "{}".format(",".join(variant.callers) or None))
            sheet.write(row, 13, "{}".format(variant.num_times_callers))
            sheet.write(row, 14, "{}".format(",".join(variant.cosmic_ids) or None))
            sheet.write(row, 15, "{}".format(variant.cosmic_data['num_samples']))
            sheet.write(row, 16, "{}".format(variant.cosmic_data['aa']))
            sheet.write(row, 17, "{}".format(variant.clinvar_data['significance']))
            sheet.write(row, 18, "{}".format(variant.clinvar_data['hgvs']))
            sheet.write(row, 19, "{}".format(variant.clinvar_data['disease']))
            sheet.write(row, 20, "{}".format(coverage_string))
            sheet.write(row, 21, "{}".format(reads_string))
            sheet.write(row, 22, "{}".format(variant.impact))
            sheet.write(row, 23, "{}".format(variant.severity))
            sheet.write(row, 24, "{}".format(variant.max_maf_all))
            sheet.write(row, 25, "{}".format(variant.min_depth))
            sheet.write(row, 26, "{}".format(variant.max_depth))
            sheet.write(row, 27, "{}".format(variant.chr))
            sheet.write(row, 28, "{}".format(variant.pos))
            sheet.write(row, 29, "{}".format(variant.end))
            sheet.write(row, 30, "{}".format(",".join(variant.rs_ids)))

            col = 31
            if 'mutect' in callers:
                sheet.write(row, col, "{}".format(variant.mutect.get('AAF') or None))
                col += 1

            if 'vardict' in callers:
                sheet.write(row, col, "{}".format(variant.vardict.get('AAF') or None))
                col += 1

            if 'freebayes' in callers:
                sheet.write(row, col, "{}".format(variant.freebayes.get('AAF') or None))
                col += 1

            if 'scalpel' in callers:
                sheet.write(row, col, "{}".format(variant.scalpel.get('AAF') or None))
                col += 1

            if 'platypus' in callers:
                sheet.write(row, col, "{}".format(variant.platypus.get('AAF') or None))
                col += 1

            if 'pindel' in callers:
                sheet.write(row, col, "{}".format(variant.pindel.get('AAF') or None))
                col += 1

            row += 1
        sheet_num += 1
    wb.save(report_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-r', '--report', help="Root name for reports (per sample)", default='report')
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    parser.add_argument('-d', '--min_depth', help='Minimum depth threshold for variant reporting', default=250.0)
    parser.add_argument('-t', '--min_somatic_var_freq', help='Minimum reportable somatic variant frequency',
                        default=0.10)
    parser.add_argument('-p', '--max_pop_freq', help='Maximum allowed population allele frequency', default=0.005)
    Job.Runner.addToilOptions(parser)
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
    else:
        auth_provider = None

    thresholds = {'min_saf': args.min_somatic_var_freq,
                  'max_maf': args.max_pop_freq,
                  'depth': args.min_depth}

    callers = ("mutect", "platypus", "vardict", "scalpel", "freebayes", "pindel")

    sys.stdout.write("Processing samples\n")
    root_job = Job.wrapJobFn(pipeline.spawn_batch_jobs, cores=1)

    for sample in samples:
        sample_job = Job.wrapJobFn(process_sample, config, sample, samples, [args.address], auth_provider,
                                   thresholds, callers, cores=1)

        root_job.addChild(sample_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)


import sys
import csv
import numpy as np
import geneimpacts

from collections import defaultdict

from variantstore import Variant
from variantstore import SampleVariant
from coveragestore import SampleCoverage
from cassandra.cqlengine import connection


def get_target_amplicons(filename):
    amplicons_list = list()
    with open(filename, "r") as bedfile:
        reader = csv.reader(bedfile, dialect='excel-tab')
        for row in reader:
            amplicons_list.append(row[3])

    return amplicons_list


def get_preferred_transcripts(transcripts_file):
    transcripts_list = list()
    with open(transcripts_file, 'r') as infile:
        reader = csv.reader(infile, dialect='excel-tab')
        for row in reader:
            if row[2] not in transcripts_list:
                transcripts_list.append(row[2])

    return transcripts_list


def get_effects(variant, annotation_keys):
    effects = []
    effects += [geneimpacts.SnpEff(e, annotation_keys) for e in variant.INFO.get("ANN").split(",")]

    return effects


def get_top_impact(effects):
    top_impact = geneimpacts.Effect.top_severity(effects)

    if isinstance(top_impact, list):
        top_impact = top_impact[0]

    return top_impact


def get_genes(effects):
    genes_list = []

    for effect in effects:
        if effect.gene not in genes_list:
            genes_list.append(effect.gene)

    return genes_list


def get_transcript_effects(effects):
    transcript_effects = dict()

    for effect in effects:
        if effect.transcript is not None:
            transcript_effects[effect.transcript] = "{gene}|{biotype}|{impact}|{exon}|{codon}|{aa}|{cons}|{so}" \
                                                    "".format(gene=effect.gene,
                                                              biotype=effect.biotype,
                                                              impact=effect.impact_severity,
                                                              exon=effect.exon,
                                                              codon=effect.codon_change,
                                                              aa=effect.aa_change,
                                                              cons=effect.top_consequence,
                                                              so=effect.so)

    return transcript_effects


def get_clinvar_info(variant, samples, sample):
    clinvar_data = dict()

    clinvar_data['significance'] = variant.INFO.get('clinvar_significance') or 'None'
    clinvar_data['pathogenic'] = variant.INFO.get('clinvar_pathogenic') or 'None'
    clinvar_data['hgvs'] = variant.INFO.get('clinvar_hgvs') or 'None'
    clinvar_data['revstatus'] = variant.INFO.get('clinvar_revstatus') or 'None'
    clinvar_data['org'] = variant.INFO.get('clinvar_org') or 'None'
    clinvar_data['disease'] = variant.INFO.get('clinvar_diseasename') or 'None'
    clinvar_data['accession'] = variant.INFO.get('clinvar_accession') or 'None'

    try:
        clinvar_data['origin'] = variant.INFO.get('clinvar_origin') or 'None'
    except IndexError:
        clinvar_data['origin'] = 'None'
        with open("{}.sample_variant_add.log".format(samples[sample]['library_name']), "a") as err:
            err.write("Problem with ClinVar origin data for variant, setting to None:\n")
            err.write("Sample: {}\t Library: {}\n".format(samples[sample]['sample_name'],
                                                          samples[sample]['library_name'], ))
            err.write("{}\n".format(variant))

    return clinvar_data


def get_cosmic_info(variant):
    cosmic_data = dict()

    cosmic_data['ids'] = variant.INFO.get('cosmic_ids') or 'None'
    cosmic_data['num_samples'] = unicode(variant.INFO.get('cosmic_numsamples')) or 'None'
    cosmic_data['cds'] = variant.INFO.get('cosmic_cds') or 'None'
    cosmic_data['aa'] = variant.INFO.get('cosmic_aa') or 'None'
    cosmic_data['gene'] = variant.INFO.get('cosmic_gene') or 'None'

    return cosmic_data


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
            'panel_amplicon': variant.INFO.get('panel_target') or "None",
            'intersect': variant.INFO.get('amplicon_intersect') or "None"}

    return data


def get_variants(config, samples, sample, library, thresholds, report_names):
    variants = SampleVariant.objects.timeout(None).filter(
        SampleVariant.reference_genome == config['genome_version'],
        SampleVariant.sample == samples[sample][library]['sample_name'],
        SampleVariant.run_id == samples[sample][library]['run_id'],
        SampleVariant.library_name == samples[sample][library]['library_name'],
        SampleVariant.max_maf_all <= thresholds['max_maf']
    ).allow_filtering()

    num_var = variants.count()
    sys.stdout.write("Retrieved {} total variants\n".format(num_var))
    with open(report_names['log'], 'a') as logfile:
        logfile.write("Retrieved {} total variants\n".format(num_var))

    ordered = variants.order_by('library_name', 'chr', 'pos', 'ref', 'alt').limit(variants.count() + 1000)

    return ordered, num_var


def get_coverage_data(target_amplicons, samples, sample, library, target_amplicon_coverage):
    reportable_amplicons = list()
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
            reportable_amplicons.append(result)
            target_amplicon_coverage[amplicon]['num_reads'] = result.num_reads
            target_amplicon_coverage[amplicon]['mean_coverage'] = result.mean_coverage

    return reportable_amplicons, target_amplicon_coverage


def setup_report_header(filename, callers):
    with open(filename, 'w') as report:
        report.write("Sample\tLibrary\tGene\tAmplicon\tRef\tAlt\tCodon\tAA\t"
                     "Max VAF\tCallers\tCOSMIC_IDs\tCOSMIC_NumSamples\tCOSMIC_AA\t"
                     "Clin_Sig\tClin_HGVS\tClin_Disease\tCoverage\tNum Reads\tImpact\tSeverity\tmax_maf_all\t"
                     "max_maf_no_fin\tmin_caller_depth\tmax_caller_depth\tChrom\tStart\tEnd\trsIDs")

        if 'mutect' in callers:
            report.write("\tMuTect_AF")

        if 'vardict' in callers:
            report.write("\tVarDict_AF")

        if 'freebayes' in callers:
            report.write("\tFreeBayes_AF")

        if 'scalpel' in callers:
            report.write("\tScalpel_AF")

        if 'platypus' in callers:
            report.write("\tPlatypus_AF")

        if 'pindel' in callers:
            report.write("\tPindel_AF")

        report.write("\n")


def classify_and_filter_variants_proj(samples, sample, library, report_names, target_amplicons, callers,
                                      ordered_variants, config, thresholds, project_variant_data, variant_count_data,
                                      gene_count_data, address, auth_provider):

    iterated = 0
    passing_variants = 0
    filtered_low_freq = 0

    tier1_pass_variants = list()
    tier1_fail_variants = list()

    vus_pass_variants = list()
    vus_fail_variants = list()

    tier4_pass_variants = list()
    tier4_fail_variants = list()

    filtered_off_target = list()
    off_target_amplicon_counts = defaultdict(int)

    category = samples[sample][library]['category']
    counted = list()

    connection.setup([address], "variantstore", auth_provider=auth_provider)
    sample_keys = samples.keys()

    for variant in ordered_variants:
        if len(variant.callers) < 2:
            continue
        if len(variant.ref) > 2 and len(variant.alt) > 2:
            continue

        iterated += 1
        variant_id = "{}:{}-{}_{}_{}_{}_{}".format(variant.chr, variant.pos, variant.end, variant.ref, variant.alt,
                                                   variant.codon_change, variant.aa_change)
        project_variant_data[variant_id]['variant'] = variant

        elements = variant.amplicon_data['amplicon'].split('_')
        gene = elements[0]

        if variant.amplicon_data['amplicon'] is 'None':
            filtered_off_target.append(variant)
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
                            cohort_count = 0.0
                            counted_samples = list()
                            for var in ordered_var:
                                if var.sample in sample_keys:
                                    if var.sample not in counted_samples:
                                        cohort_count += 1
                                        counted_samples.append(var.sample)
                            fraction = cohort_count / len(sample_keys)
                            if fraction < 0.5:
                                counted.append(variant_id)
                                project_variant_data[variant_id][category] += 1
                                variant_count_data[sample]['pass_count'] += 1
                                gene_count_data[sample][gene] += 1
                                if variant.ref == 'C' and variant.alt == 'T':
                                    variant_count_data[sample]['CT_count'] += 1

                                if variant.severity == 'HIGH' or 'pathogenic' in variant.clinvar_data['pathogenic']:
                                    variant_count_data[sample]['high_impact_pathogenic'] += 1
                                elif variant.severity == 'MED':
                                    variant_count_data[sample]['med'] += 1
                                elif variant.severity == 'LOW':
                                    variant_count_data[sample]['low'] += 1
                                else:
                                    sys.stderr.write("ERROR: Cannot classify variant {}\n".format(variant_id))

                                # Putting in to Tier1 based on COSMIC
                                if variant.cosmic_ids:
                                    tier1_pass_variants.append(variant)
                                    project_variant_data[variant_id]['tier1_pass'] += 1
                                    passing_variants += 1
                                    continue

                                # Putting in to Tier1 based on ClinVar not being None or Benign
                                if variant.clinvar_data['pathogenic'] != 'None':
                                    if variant.clinvar_data['pathogenic'] != 'benign':
                                        if variant.clinvar_data['pathogenic'] != 'likely-benign':
                                            tier1_pass_variants.append(variant)
                                            project_variant_data[variant_id]['tier1_pass'] += 1
                                            passing_variants += 1
                                            continue

                                if variant.severity == 'MED' or variant.severity == 'HIGH':
                                    vus_pass_variants.append(variant)
                                    project_variant_data[variant_id]['vus_pass'] += 1
                                    passing_variants += 1
                                    continue
                                else:
                                    tier4_pass_variants.append(variant)
                                    project_variant_data[variant_id]['tier4_pass'] += 1
                                    passing_variants += 1
                                    continue
                        else:
                            # sys.stderr.write("WARNING: Duplicate variant, skipping: {}\n".format(variant_id))
                            with open("{}_Duplicates.log".format(sample), 'a') as duplicates:
                                duplicates.write("{}\n".format(variant_id))
            else:
                filtered_off_target.append(variant)
                off_target_amplicon_counts[variant.amplicon_data['amplicon']] += 1

    sys.stdout.write("Iterated through {} variants\n".format(iterated))
    with open(report_names['log'], 'a') as logfile:
        logfile.write("Iterated through {} variants\n".format(iterated))
        logfile.write("---------------------------------------------\n")
        logfile.write("{}\n".format(library))
        logfile.write(
            "Sending {} variants to reporting (filtered {} off-target  and {} low frequency/low depth variants)"
            "\n".format(passing_variants, filtered_low_freq, len(filtered_off_target)))
        logfile.write("---------------------------------------------\n")
        logfile.write("Off Target Amplicon\tCounts\n")

        sys.stdout.write("---------------------------------------------\n")
        sys.stdout.write("{}\n".format(library))
        sys.stdout.write(
            "Sending {} variants to reporting (filtered {} off-target  and {} low frequency/low depth variants)"
            "\n".format(passing_variants, filtered_low_freq, len(filtered_off_target)))
        sys.stdout.write("---------------------------------------------\n")
        sys.stdout.write("Off Target Amplicon\tCounts\n")

        for off_target in off_target_amplicon_counts:
            logfile.write("{}\t{}\n".format(off_target, off_target_amplicon_counts[off_target]))
            # sys.stdout.write("{}\t{}\n".format(off_target, off_target_amplicon_counts[off_target]))

    return tier1_pass_variants, tier1_fail_variants, vus_pass_variants, vus_fail_variants, tier4_pass_variants, \
           tier4_fail_variants, filtered_off_target, project_variant_data, variant_count_data, gene_count_data


def write_reports(report_names, samples, sample, library, filtered_var_data, ordered_variants,
                  target_amplicon_coverage, reportable_amplicons, num_var, thresholds, callers):

    tier1_pass_variants, tier1_fail_variants, vus_pass_variants, vus_fail_variants, tier4_pass_variants, \
    tier4_fail_variants, filtered_off_target, project_variant_data, variant_count_data, \
    gene_count_data = filtered_var_data

    with open(report_names['coverage'], "a") as coverage_report:
        coverage_report.write("Library:\t{}\n".format(samples[sample][library]['library_name']))
        coverage_report.write("Run ID:\t{}\n".format(samples[sample][library]['run_id']))
        coverage_report.write("Report:\t{}\n".format(samples[sample][library]['report']))
        coverage_report.write("Min Somatic Allele Frequency:\t{}\n".format(thresholds['min_saf']))
        coverage_report.write("Max Population Germline Allele Frequency:\t{}\n".format(thresholds['max_maf']))
        coverage_report.write("---------------------------------------------\n")
        coverage_report.write("Number query return variants\t{}\n".format(num_var))
        coverage_report.write("Tier1 Pass Variants\t{}\n".format(len(tier1_pass_variants)))
        coverage_report.write("Tier1 Fail Variants\t{}\n".format(len(tier1_fail_variants)))

        coverage_report.write("VUS Pass Variants\t{}\n".format(len(vus_pass_variants)))
        coverage_report.write("VUS Fail Variants\t{}\n".format(len(vus_fail_variants)))

        coverage_report.write("Tier4 Pass Variants\t{}\n".format(len(tier4_pass_variants)))
        coverage_report.write("Tier4 Fail Variants\t{}\n".format(len(tier4_fail_variants)))

        coverage_report.write("Off Target/Non-Reportable Amplicon Variants\t{}\n".format(len(filtered_off_target)))
        coverage_report.write("---------------------------------------------\n")
        coverage_report.write("Sample\tLibrary\tAmplicon\tNum Reads\tCoverage\n")
        for amplicon in reportable_amplicons:
            coverage_report.write("{}\t{}\t{}\t{}\t{}\n".format(amplicon.sample,
                                                                amplicon.library_name,
                                                                amplicon.amplicon,
                                                                amplicon.num_reads,
                                                                amplicon.mean_coverage))

        write_report(report_names['tier1_pass'], tier1_pass_variants, target_amplicon_coverage, callers)
        write_report(report_names['tier1_fail'], tier1_fail_variants, target_amplicon_coverage, callers)

        write_report(report_names['vus_pass'], vus_pass_variants, target_amplicon_coverage, callers)
        write_report(report_names['vus_fail'], vus_fail_variants, target_amplicon_coverage, callers)

        write_report(report_names['tier4_pass'], tier4_pass_variants, target_amplicon_coverage, callers)
        write_report(report_names['tier4_fail'], tier4_fail_variants, target_amplicon_coverage, callers)

        write_report(report_names['all_ordered'], tier4_fail_variants, target_amplicon_coverage, callers)


def write_report(filename, variants, target_amplicon_coverage, callers):
    with open(filename, 'a') as report:
        for variant in variants:
            try:
                report.write("{sample}\t{library}\t{gene}\t{amp}\t{ref}\t{alt}\t{codon}\t{aa}\t"
                             "{max_som_aaf}\t{callers}\t{cosmic}\t{cosmic_nsamples}\t{cosmic_aa}\t"
                             "{csig}\t{hgvs}\t{cdis}\t{cov}\t{reads}\t{impact}\t{severity}\t{max_maf_all}\t"
                             "{max_maf_no_fin}\t{min_depth}\t{max_depth}\t{chr}\t{start}\t{end}\t{rsids}"
                             "".format(sample=variant.sample,
                                       library=variant.library_name,
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
                                       severity=variant.severity,
                                       max_maf_all=variant.max_maf_all,
                                       max_maf_no_fin=variant.max_maf_no_fin,
                                       max_som_aaf=variant.max_som_aaf,
                                       min_depth=variant.min_depth,
                                       max_depth=variant.max_depth,
                                       callers=",".join(variant.callers) or None))
            except KeyError:
                sys.stderr.write("Could not write variant to report, KeyError with missing data\n")
                print variant
                continue

            if 'mutect' in callers:
                report.write("\t{maf}"
                             "".format(maf=variant.mutect.get('AAF') or None))

            if 'vardict' in callers:
                report.write("\t{vaf}"
                             "".format(vaf=variant.vardict.get('AAF') or None))

            if 'freebayes' in callers:
                report.write("\t{faf}"
                             "".format(faf=variant.freebayes.get('AAF') or None))

            if 'scalpel' in callers:
                report.write("\t{saf}"
                             "".format(saf=variant.scalpel.get('AAF') or None))

            if 'platypus' in callers:
                report.write("\t{plaf}"
                             "".format(plaf=variant.platypus.get('AAF') or None))

            if 'pindel' in callers:
                report.write("\t{paf}"
                             "".format(paf=variant.pindel.get('AAF') or None))

            report.write("\n")


def write_sample_variant_report(report_root, sample, variants, target_amplicon_coverage, callers):
    with open("{}_variant_{}.txt".format(sample, report_root), 'w') as report:
        report.write("Sample\tLibrary\tGene\tAmplicon\tRef\tAlt\tCodon\tAA\t"
                     "max_somatic_aaf\tCallers\tCOSMIC_IDs\tCOSMIC_NumSamples\tCOSMIC_AA\t"
                     "Clin_Sig\tClin_HGVS\tClin_Disease\t"
                     "Coverage\tNum Reads\tImpact\tSeverity\tmax_maf_all\tmax_maf_no_fin\t"
                     "min_caller_depth\tmax_caller_depth\tChrom\tStart\tEnd\trsIDs")

        if 'mutect' in callers:
            report.write("\tMuTect_AF")

        if 'vardict' in callers:
            report.write("\tVarDict_AF")

        if 'freebayes' in callers:
            report.write("\tFreeBayes_AF")

        if 'scalpel' in callers:
            report.write("\tScalpel_AF")

        if 'platypus' in callers:
            report.write("\tPlatypus_AF")

        if 'pindel' in callers:
            report.write("\tPindel_AF")

        report.write("\n")

        num_reported = 0

        for variant in variants:
            # print variant
            if any(caller in ("freebayes", "pindel") for caller in variant.callers):
                # sys.stderr.write("Freebayes or Pindel detected in callers: {}\n".format(variant.callers))
                if any(caller in ("mutect", "scalpel", "vardict", "platypus") for caller in variant.callers):
                    # Other caller also called this variant
                    pass
                else:
                    # Freebayes or Pindel called only
                    if variant.cosmic_ids:
                        # Pindel/FreeBayes only but COSMIC IDs
                        pass
                    elif variant.clinvar_data['pathogenic'] != 'None':
                        # Pindel/FreeBayes only but Clinvar data
                        pass
                    else:
                        # Freebayes or Pindel only, no cosmic or clinvar data
                        continue
            try:
                report.write("{sample}\t{library}\t{gene}\t{amp}\t{ref}\t{alt}\t{codon}\t{aa}\t"
                             "{max_som_aaf}\t{callers}\t{cosmic}\t{cosmic_nsamples}\t{cosmic_aa}\t{csig}\t{hgvs}\t"
                             "{cdis}\t{cov}\t{reads}\t{impact}\t{severity}\t{max_maf_all}\t{max_maf_no_fin}\t"
                             "{min_depth}\t{max_depth}\t{chr}\t{start}\t{end}\t{rsids}"
                             "".format(sample=variant.sample,
                                       library=variant.library_name,
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
                                       severity=variant.severity,
                                       max_maf_all=variant.max_maf_all,
                                       max_maf_no_fin=variant.max_maf_no_fin,
                                       max_som_aaf=variant.max_som_aaf,
                                       min_depth=variant.min_depth,
                                       max_depth=variant.max_depth,
                                       callers=",".join(variant.callers) or None))
            except KeyError:
                sys.stderr.write("Could not write variant to report, KeyError with missing data\n")
                print variant
                continue

            num_reported += 1

            if 'mutect' in callers:
                report.write("\t{maf}"
                             "".format(maf=variant.mutect.get('AAF') or None))

            if 'vardict' in callers:
                report.write("\t{vaf}"
                             "".format(vaf=variant.vardict.get('AAF') or None))

            if 'freebayes' in callers:
                report.write("\t{faf}"
                             "".format(faf=variant.freebayes.get('AAF') or None))

            if 'scalpel' in callers:
                report.write("\t{saf}"
                             "".format(saf=variant.scalpel.get('AAF') or None))

            if 'platypus' in callers:
                report.write("\t{plaf}"
                             "".format(plaf=variant.platypus.get('AAF') or None))

            if 'pindel' in callers:
                report.write("\t{paf}"
                             "".format(paf=variant.pindel.get('AAF') or None))

            report.write("\n")

        sys.stdout.write("Wrote {} variants in report\n".format(num_reported))


def write_sample_variant_report_no_caller_filter(report_root, sample, variants, target_amplicon_coverage, callers):
    with open("{}.{}.txt".format(sample, report_root), 'w') as report:
        report.write("Sample\tLibrary\tGene\tAmplicon\tRef\tAlt\tCodon\tAA\t"
                     "max_somatic_aaf\tCallers\tCOSMIC_IDs\tCOSMIC_NumSamples\tCOSMIC_AA\t"
                     "Clin_Sig\tClin_HGVS\tClin_Disease\t"
                     "Coverage\tNum Reads\tImpact\tSeverity\tmax_maf_all\tmax_maf_no_fin\t"
                     "min_caller_depth\tmax_caller_depth\tChrom\tStart\tEnd\trsIDs")

        if 'mutect' in callers:
            report.write("\tMuTect_AF")

        if 'vardict' in callers:
            report.write("\tVarDict_AF")

        if 'freebayes' in callers:
            report.write("\tFreeBayes_AF")

        if 'scalpel' in callers:
            report.write("\tScalpel_AF")

        if 'platypus' in callers:
            report.write("\tPlatypus_AF")

        if 'pindel' in callers:
            report.write("\tPindel_AF")

        report.write("\n")

        num_reported = 0

        for variant in variants:
            report.write("{sample}\t{library}\t{gene}\t{amp}\t{ref}\t{alt}\t{codon}\t{aa}\t"
                         "{max_som_aaf}\t{callers}\t{cosmic}\t{cosmic_nsamples}\t{cosmic_aa}\t{csig}\t{hgvs}\t{cdis}\t"
                         "{cov}\t{reads}\t{impact}\t{severity}\t{max_maf_all}\t{max_maf_no_fin}\t"
                         "{min_depth}\t{max_depth}\t{chr}\t{start}\t{end}\t{rsids}"
                         "".format(sample=variant.sample,
                                   library=variant.library_name,
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
                                   severity=variant.severity,
                                   max_maf_all=variant.max_maf_all,
                                   max_maf_no_fin=variant.max_maf_no_fin,
                                   max_som_aaf=variant.max_som_aaf,
                                   min_depth=variant.min_depth,
                                   max_depth=variant.max_depth,
                                   callers=",".join(variant.callers) or None))
            num_reported += 1

            if 'mutect' in callers:
                report.write("\t{maf}"
                             "".format(maf=variant.mutect.get('AAF') or None))

            if 'vardict' in callers:
                report.write("\t{vaf}"
                             "".format(vaf=variant.vardict.get('AAF') or None))

            if 'freebayes' in callers:
                report.write("\t{faf}"
                             "".format(faf=variant.freebayes.get('AAF') or None))

            if 'scalpel' in callers:
                report.write("\t{saf}"
                             "".format(saf=variant.scalpel.get('AAF') or None))

            if 'platypus' in callers:
                report.write("\t{plaf}"
                             "".format(plaf=variant.platypus.get('AAF') or None))

            if 'pindel' in callers:
                report.write("\t{paf}"
                             "".format(paf=variant.pindel.get('AAF') or None))

            report.write("\n")

        sys.stdout.write("Wrote {} variants in report\n".format(num_reported))

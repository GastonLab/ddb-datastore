import sys
import csv
import geneimpacts

from collections import defaultdict
from variantstore import SampleVariant
from coveragestore import SampleCoverage


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


def get_variants(config, samples, sample, library, thresholds, report_root):
    variants = SampleVariant.objects.timeout(None).filter(
        SampleVariant.reference_genome == config['genome_version'],
        SampleVariant.sample == samples[sample][library]['sample_name'],
        SampleVariant.run_id == samples[sample][library]['run_id'],
        SampleVariant.library_name == samples[sample][library]['library_name'],
        SampleVariant.max_maf_all <= thresholds['max_maf'],
    ).allow_filtering()

    sys.stdout.write("Retrieved {} total variants\n".format(variants.count()))
    with open("{}.{}.log".format(sample, report_root), 'a') as logfile:
        logfile.write("Retrieved {} total variants\n".format(variants.count()))

    ordered = variants.order_by('library_name', 'chr', 'pos', 'ref', 'alt').limit(variants.count() + 1000)

    return ordered


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


def variant_filter(variant, thresholds):
    flag = False
    info = dict()

    # The ClinVar thing needs some work to be processed appropriately
    if variant.in_cosmic or variant.clinvar_data['significance'] is not 'benign':
        flag = True

    if variant.clinvar_data['significance'] is 'None':
        info['clinvar'] = "Unknown"
    elif variant.clinvar_data['significance'] is not 'benign':
        info['clinvar'] = "Not Benign"
    else:
        info['clinvar'] = "Potentially Benign"

    if variant.max_maf_all < thresholds['max_maf']:
        info['max_maf'] = "Not Common"
    else:
        info['max_maf'] = "Common"

    if variant.amplicon_data['amplicon'] is not 'None':
        info['dual'] = True
    else:
        info['dual'] = False

    return flag, info


def filter_variants_on_target(ordered_variants, target_amplicons):
    on_target_variants = list()
    filtered_no_amplicon = list()
    filtered_non_target_amplicon = list()
    off_target_amplicons = defaultdict(int)

    for variant in ordered_variants:
        if variant.amplicon_data['amplicon']:
            amplicons = variant.amplicon_data['amplicon'].split(',')
            for amplicon in amplicons:
                if amplicon in target_amplicons:
                    on_target_variants.append(variant)
                    break
                else:
                    filtered_non_target_amplicon.append(variant)
                    off_target_amplicons[amplicon] += 1
        elif variant.amplicon_data['amplicon'] is 'None':
            filtered_no_amplicon.append(variant)
            off_target_amplicons[amplicon] += 1
        else:
            filtered_no_amplicon.append(variant)
            off_target_amplicons[amplicon] += 1

    return on_target_variants, filtered_no_amplicon, filtered_non_target_amplicon, off_target_amplicons


def setup_report_header(filename, callers):
    with open(filename, 'w') as report:
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


def classify_and_filter_variants(sample, library, report_root, target_amplicons, callers, ordered_variants, thresholds):

    iterated = 0
    passing_variants = 0
    filtered_low_freq = 0

    tier1_pass_variants = list()
    tier1_fail_variants = list()

    tier2_pass_variants = list()
    tier2_fail_variants = list()

    vus_pass_variants = list()
    vus_fail_variants = list()

    tier4_pass_variants = list()
    tier4_fail_variants = list()

    filtered_off_target = list()
    off_target_amplicon_counts = defaultdict(int)

    for variant in ordered_variants:
        iterated += 1
        if variant.amplicon_data['amplicon']:
            amplicons = variant.amplicon_data['amplicon'].split(',')
            for amplicon in amplicons:
                if amplicon in target_amplicons:
                    for caller in callers:
                        if caller in variant.callers:

                            # Putting in to Tier1 based on COSMIC
                            if variant.cosmic_ids:
                                if variant.max_som_aaf < thresholds['min_saf']:
                                    tier1_fail_variants.append(variant)
                                    filtered_low_freq += 1
                                else:
                                    tier1_pass_variants.append(variant)
                                    passing_variants += 1
                                break

                            # Putting in to Tier1 based on ClinVar not being None or Benign
                            if variant.clinvar_data['pathogenic'] != 'None':
                                if variant.clinvar_data['pathogenic'] != 'benign':
                                    if variant.clinvar_data['pathogenic'] != 'likely-benign':
                                        if variant.max_som_aaf < thresholds['min_saf']:
                                            tier1_fail_variants.append(variant)
                                            filtered_low_freq += 1
                                        else:
                                            tier1_pass_variants.append(variant)
                                            passing_variants += 1
                                break

                            if variant.severity == 'MED' or variant.severity == 'HIGH':
                                if variant.max_som_aaf < thresholds['min_saf']:
                                    vus_fail_variants.append(variant)
                                    filtered_low_freq += 1
                                else:
                                    vus_pass_variants.append(variant)
                                    passing_variants += 1
                            else:
                                if variant.max_som_aaf < thresholds['min_saf']:
                                    tier4_fail_variants.append(variant)
                                    filtered_low_freq += 1
                                else:
                                    tier4_pass_variants.append(variant)
                                    passing_variants += 1
                else:
                    filtered_off_target.append(variant)
                    off_target_amplicon_counts[amplicon] += 1
        elif variant.amplicon_data['amplicon'] is 'None':
            filtered_off_target.append(variant)
            off_target_amplicon_counts[amplicon] += 1
        else:
            filtered_off_target.append(variant)
            off_target_amplicon_counts[amplicon] += 1

    sys.stdout.write("Iterated through {} variants\n".format(iterated))
    with open("{}.{}.log".format(sample, report_root), 'a') as logfile:
        logfile.write("Iterated through {} variants\n".format(iterated))

    with open("{}.{}.log".format(sample, report_root), 'a') as logfile:
        logfile.write("---------------------------------------------\n")
        logfile.write("{}\n".format(library))
        logfile.write(
            "Sending {} variants to reporting (filtered {} off-target  and {} low-frequency variants)"
            "\n".format(passing_variants, filtered_low_freq, len(filtered_off_target)))
        logfile.write("---------------------------------------------\n")
        logfile.write("Off Target Amplicon\tCounts\n")

        sys.stdout.write("---------------------------------------------\n")
        sys.stdout.write("{}\n".format(library))
        sys.stdout.write(
            "Sending {} variants to reporting (filtered {} off-target  and {} low-frequency variants)"
            "\n".format(passing_variants, filtered_low_freq, len(filtered_off_target)))
        sys.stdout.write("---------------------------------------------\n")
        sys.stdout.write("Off Target Amplicon\tCounts\n")

        for off_target in off_target_amplicon_counts:
            logfile.write("{}\t{}\n".format(off_target, off_target_amplicon_counts[off_target]))
            sys.stdout.write("{}\t{}\n".format(off_target, off_target_amplicon_counts[off_target]))

    return tier1_pass_variants, tier1_fail_variants, vus_pass_variants, vus_fail_variants, tier4_pass_variants, \
           tier4_fail_variants, filtered_off_target, off_target_amplicon_counts


def write_reports(report_root, samples, sample, library, filtered_var_data, target_amplicon_coverage,
                  reportable_amplicons, thresholds, callers):
    tier1_pass_variants, tier1_fail_variants, vus_pass_variants, vus_fail_variants, tier4_pass_variants, \
    tier4_fail_variants, filtered_off_target, off_target_amplicon_counts = filtered_var_data

    with open("{}_coverage_{}.txt".format(sample, report_root), "a") as coverage_report:
        coverage_report.write("Library:\t{}\n".format(samples[sample][library]['library_name']))
        coverage_report.write("Run ID:\t{}\n".format(samples[sample][library]['run_id']))
        coverage_report.write("Report:\t{}\n".format(samples[sample][library]['report']))
        coverage_report.write("Min Somatic Allele Frequency:\t{}\n".format(thresholds['min_saf']))
        coverage_report.write("Max Population Germline Allele Frequency:\t{}\n".format(thresholds['max_maf']))
        coverage_report.write("---------------------------------------------\n")
        coverage_report.write("Sample\tLibrary\tAmplicon\tNum Reads\tCoverage\n")
        for amplicon in reportable_amplicons:
            coverage_report.write("{}\t{}\t{}\t{}\t{}\n".format(amplicon.sample,
                                                                amplicon.library_name,
                                                                amplicon.amplicon,
                                                                amplicon.num_reads,
                                                                amplicon.mean_coverage))

        write_report("{}_tier1_pass_variants_{}.txt".format(sample, report_root), tier1_pass_variants,
                     target_amplicon_coverage, callers)
        write_report("{}_tier1_fail_variants_{}.txt".format(sample, report_root), tier1_fail_variants,
                     target_amplicon_coverage, callers)

        write_report("{}_vus_pass_variants_{}.txt".format(sample, report_root), vus_pass_variants,
                     target_amplicon_coverage, callers)
        write_report("{}_vus_fail_variants_{}.txt".format(sample, report_root), vus_fail_variants,
                     target_amplicon_coverage, callers)

        write_report("{}_tier4_pass_variants_{}.txt".format(sample, report_root), tier4_pass_variants,
                     target_amplicon_coverage, callers)
        write_report("{}_tier4_fail_variants_{}.txt".format(sample, report_root), tier4_fail_variants,
                     target_amplicon_coverage, callers)


def write_report(filename, variants, target_amplicon_coverage, callers):
    with open(filename, 'a') as report:
        for variant in variants:
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

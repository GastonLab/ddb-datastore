import sys
import geneimpacts
from pybedtools import BedTool


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
            transcript_effects[effect.transcript] = "{biotype}|{effect}".format(biotype=effect.biotype,
                                                                                effect=effect.impact_severity)

    return transcript_effects


def get_clinvar_info(variant):
    clinvar_data = dict()

    clinvar_data['significance'] = variant.INFO.get('clinvar_significance') or 'None'
    clinvar_data['pathogenic'] = variant.INFO.get('clinvar_pathogenic') or 'None'
    clinvar_data['hgvs'] = variant.INFO.get('clinvar_hgvs') or 'None'
    clinvar_data['revstatus'] = variant.INFO.get('clinvar_revstatus') or 'None'
    clinvar_data['origin'] = variant.INFO.get('clinvar_origin') or 'None'
    clinvar_data['org'] = variant.INFO.get('clinvar_org') or 'None'
    clinvar_data['disease'] = variant.INFO.get('clinvar_diseasename') or 'None'
    clinvar_data['accession'] = variant.INFO.get('clinvar_accession') or 'None'

    return clinvar_data


def get_cosmic_info(variant):
    cosmic_data = dict()

    cosmic_data['ids'] = variant.INFO.get('cosmic_ids') or 'None'
    cosmic_data['num_samples'] = variant.INFO.get('cosmic_numsamples') or 'None'
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


def write_sample_variant_report(report_root, sample, variants, target_amplicon_coverage, callers):
    with open("{}.{}.txt".format(sample, report_root), 'w') as report:
        report.write("Chrom\tStart\tEnd\tGene\tRef\tAlt\tCodon\tAA\trsIDs\tAmplicon\t"
                     "COSMIC_IDs\tCOSMIC_NumSamples\tClin_Sig\tClin_HGVS\tClin_Disease\t"
                     "Biotype\tImpact\tImpact SO\tSeverity\tmax_maf_all\tmax_maf_no_fin\tmax_somatic_aaf\t"
                     "NumReads\tCoverage\tmin_caller_depth\tmax_caller_depth\tCallers")

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
                    # sys.stderr.write("Other caller also called this variant\n")
                    pass
                else:
                    # sys.stderr.write("Freebayes or Pindel called only\n")
                    if variant.cosmic_ids:
                        sys.stderr.write("Pindel/FreeBayes only but COSMIC IDs ({}) data\n".format(variant.cosmic_ids))
                    elif variant.clinvar_data['significance'] is not None:
                        sys.stderr.write("Pindel/FreeBayes only but Clinvar data ({})\n".format(variant.clinvar_data))
                    else:
                        sys.stderr.write("Freebayes or Pindel only, no cosmic or clinvar data. Skipping...\n")
                        continue
            report.write("{chr}\t{start}\t{end}\t{gene}\t{ref}\t{alt}\t{codon}\t{aa}\t{rsids}\t"
                         "{amp}\t{cosmic}\t{cosmic_nsamples}\t{csig}\t{hgvs}\t{cdis}\t{biotype}\t"
                         "{impact}\t{impact_so}\t{severity}\t{max_maf_all}\t{max_maf_no_fin}\t{max_som_aaf}\t"
                         "{num_reads}\t{coverage}\t{min_depth}\t{max_depth}\t{callers}"
                         "".format(chr=variant.chr, start=variant.pos, end=variant.end,
                                   gene=variant.gene, ref=variant.ref, alt=variant.alt,
                                   codon=variant.codon_change, aa=variant.aa_change,
                                   rsids=",".join(variant.rs_ids),
                                   cosmic=",".join(variant.cosmic_ids) or None,
                                   cosmic_nsamples=variant.cosmic_data['num_samples'],
                                   amp=variant.amplicon_data['amplicon'],
                                   csig=variant.clinvar_data['significance'],
                                   hgvs=variant.clinvar_data['hgvs'],
                                   cdis=variant.clinvar_data['disease'],
                                   biotype=variant.biotype, impact=variant.impact,
                                   impact_so=variant.impact_so, severity=variant.severity,
                                   max_maf_all=variant.max_maf_all,
                                   max_maf_no_fin=variant.max_maf_no_fin,
                                   max_som_aaf=variant.max_som_aaf,
                                   num_reads=target_amplicon_coverage[variant.amplicon_data['amplicon']][
                                       'num_reads'] or None,
                                   coverage=target_amplicon_coverage[variant.amplicon_data['amplicon']][
                                       'mean_coverage'] or None,
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

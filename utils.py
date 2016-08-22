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


def variant_filter(variant, callers, thresholds):
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


def write_sample_variant_report(report_root, sample, variants, callers, thresholds):
    with open("{}.{}.txt".format(sample, report_root), 'w') as report:
        report.write("Chrom\tStart\tEnd\tGene\tRef\tAlt\tExon\tCodon\tAA\trsIDs\tAmplicon\t"
                     "COSMIC_IDs\tCOSMIC_NumSamples\tCOSMIC_AA\tClin_Sig\tClin_Pathogenic\tClin_HGVS\tClin_Disease\t"
                     "Biotype\tImpact\tImpact SO\tSeverity\tmax_maf_all\tmax_maf_no_fin\tmax_somatic_aaf\tmin_depth\t"
                     "max_depth\tCallers\t")

        if 'mutect' in callers:
            report.write("MuTect_FILTER\tMuTect_Multiallelic\tMuTect_DP\tMuTect_AD\tMuTect_AF\t")

        if 'vardict' in callers:
            report.write("VarDict_FILTER\tVarDict_Multiallelic\tVarDict_DP\tVarDict_AD\tVarDict_AF\t")

        if 'freebayes' in callers:
            report.write("FreeBayes_FILTER\tFreeBayes_Multiallelic\tFreeBayes_DP\tFreeBayes_AF\tFreeBayes_RO\t"
                         "FreeBayes_AO\t")

        if 'scalpel' in callers:
            report.write("Scalpel_FILTER\tScalpel_Multiallelic\tScalpel_DP\tScalpel_AD\tScalpel_AF\t")

        if 'platypus' in callers:
            report.write("Platypus_FILTER\tPlatypus_Multiallelic\tPlatypus_DP\tPlatypus_AD\tPlatypus_AF\t")

        if 'pindel' in callers:
            report.write("Pindel_FILTER\tPindel_Multiallelic\tPindel_DP\tPindel_AD\tPindel_AF")

        report.write("\n")

        for variant in variants:
            report.write("{chr}\t{start}\t{end}\t{gene}\t{ref}\t{alt}\t{exon}\t{codon}\t{aa}\t{rsids}\t"
                         "{amp}\t{cosmic}\t{cosmic_nsamples}\t{cosmic_aa}\t"
                         "{csig}\t{cpath}\t{hgvs}\t{cdis}\t{biotype}\t{impact}\t{impact_so}\t{severity}\t"
                         "{max_maf_all}\t{max_maf_no_fin}\t{max_som_aaf}\t{min_depth}\t{max_depth}\t{callers}"
                         "".format(sample=variant.sample, library=variant.library_name, run_id=variant.run_id,
                                   chr=variant.chr, start=variant.pos, end=variant.end,
                                   gene=variant.gene, ref=variant.ref, alt=variant.alt, exon=variant.exon,
                                   codon=variant.codon_change, aa=variant.aa_change,
                                   rsids=",".join(variant.rs_ids),
                                   cosmic=",".join(variant.cosmic_ids) or None,
                                   cosmic_nsamples=variant.cosmic_data['num_samples'],
                                   cosmic_aa=variant.cosmic_data['aa'],
                                   amp=variant.amplicon_data['amplicon'],
                                   csig=variant.clinvar_data['significance'],
                                   cpath=variant.clinvar_data['pathogenic'],
                                   hgvs=variant.clinvar_data['hgvs'],
                                   cdis=variant.clinvar_data['disease'],
                                   biotype=variant.biotype, impact=variant.impact,
                                   impact_so=variant.impact_so, severity=variant.severity,
                                   max_maf_all=variant.max_maf_all,
                                   max_maf_no_fin=variant.max_maf_no_fin,
                                   max_som_aaf=variant.max_som_aaf,
                                   min_depth=variant.min_depth,
                                   max_depth=variant.max_depth,
                                   callers=",".join(variant.callers) or None))

            if 'mutect' in callers:
                report.write("\t{mfilter}\t{mmulti}\t{mdp}\t{mad}\t{maf}"
                             "".format(mfilter=variant.mutect.get('FILTER') or None,
                                       mmulti=variant.mutect.get('MULTIALLELIC') or None,
                                       mdp=variant.mutect.get('GTF_DP') or None,
                                       mad=variant.mutect.get('GTF_AD') or None,
                                       maf=variant.mutect.get('AAF') or None))

            if 'vardict' in callers:
                report.write("\t{vfilter}\t{vmulti}\t{vdp}\t{vad}\t{vaf}"
                             "".format(vfilter=variant.vardict.get('FILTER') or None,
                                       vmulti=variant.vardict.get('MULTIALLELIC') or None,
                                       vdp=variant.vardict.get('DP') or None,
                                       vad=variant.vardict.get('AD') or None,
                                       vaf=variant.vardict.get('AAF') or None))

            if 'freebayes' in callers:
                report.write("\t{ffilter}\t{fmulti}\t{fdp}\t{faf}\t{fro}\t{fao}"
                             "".format(ffilter=variant.freebayes.get('FILTER') or None,
                                       fmulti=variant.freebayes.get('MULTIALLELIC') or None,
                                       fdp=variant.freebayes.get('DP') or None,
                                       faf=variant.freebayes.get('AAF') or None,
                                       fro=variant.freebayes.get('RO') or None,
                                       fao=variant.freebayes.get('AO') or None))

            if 'scalpel' in callers:
                report.write("\t{sfilter}\t{smulti}\t{sdp}\t{sad}\t{saf}"
                             "".format(sfilter=variant.scalpel.get('FILTER') or None,
                                       smulti=variant.scalpel.get('MULTIALLELIC') or None,
                                       sdp=variant.scalpel.get('GTF_DP') or None,
                                       sad=variant.scalpel.get('GTF_AD') or None,
                                       saf=variant.scalpel.get('AAF') or None))

            if 'platypus' in callers:
                report.write("\t{plfilter}\t{plmulti}\t{pldp}\t{plad}\t{plaf}"
                             "".format(plfilter=variant.platypus.get('FILTER') or None,
                                       plmulti=variant.platypus.get('MULTIALLELIC') or None,
                                       plad=variant.platypus.get('TR') or None,
                                       pldp=variant.platypus.get('TC') or None,
                                       plaf=variant.platypus.get('AAF') or None))

            if 'pindel' in callers:
                report.write("\t{pfilter}\t{pmulti}\t{pdp}\t{pad}\t{paf}"
                             "".format(pfilter=variant.pindel.get('FILTER') or None,
                                       pmulti=variant.pindel.get('MULTIALLELIC') or None,
                                       pdp=variant.pindel.get('DP') or None,
                                       pad=variant.pindel.get('GTF_AD') or None,
                                       paf=variant.pindel.get('AAF') or None))

            report.write("\n")


def write_variant_report(report_root, variants, callers):
    with open("{}.txt".format(report_root), 'w') as report:
        report.write("Sample\tLibrary\tRunID\tChrom\tStart\tEnd\tGene\tRef\tAlt\tExon\tCodon\tAA\trsIDs\t"
                     "Amplicon\tCOSMIC_IDs\tCOSMIC_NumSamples\tCOSMIC_AA\tClin_Sig\tClin_Pathogenic\tClin_HGVS\t"
                     "Clin_Disease\tBiotype\tImpact\tImpact SO\tSeverity\tmax_maf_all\tmax_maf_no_fin\t"
                     "max_somatic_aaf\tmin_depth\tmax_depth\tCallers\t")

        if 'mutect' in callers:
            report.write("MuTect_FILTER\tMuTect_Multiallelic\tMuTect_DP\tMuTect_AD\tMuTect_AF\t")

        if 'vardict' in callers:
            report.write("VarDict_FILTER\tVarDict_Multiallelic\tVarDict_DP\tVarDict_AD\tVarDict_AF\t")

        if 'freebayes' in callers:
            report.write("FreeBayes_FILTER\tFreeBayes_Multiallelic\tFreeBayes_DP\tFreeBayes_AF\tFreeBayes_RO\t"
                         "FreeBayes_AO\t")

        if 'scalpel' in callers:
            report.write("Scalpel_FILTER\tScalpel_Multiallelic\tScalpel_DP\tScalpel_AD\tScalpel_AF\t")

        if 'platypus' in callers:
            report.write("Platypus_FILTER\tPlatypus_Multiallelic\tPlatypus_DP\tPlatypus_AD\tPlatypus_AF\t")

        if 'pindel' in callers:
            report.write("Pindel_FILTER\tPindel_Multiallelic\tPindel_DP\tPindel_AD\tPindel_AF")

        report.write("\n")

        for variant in variants:
            report.write("{sample}\t{library}\t{run_id}\t{chr}\t{start}\t{end}\t{gene}\t{ref}\t{alt}\t{exon}\t"
                         "{codon}\t{aa}\t{rsids}\t{amp}\t{cosmic}\t{cosmic_nsamples}\t{cosmic_aa}\t"
                         "{csig}\t{cpath}\t{hgvs}\t{cdis}\t{biotype}\t{impact}\t{impact_so}\t{severity}\t"
                         "{max_maf_all}\t{max_maf_no_fin}\t{max_som_aaf}\t{min_depth}\t{max_depth}\t{callers}"
                         "".format(sample=variant.sample, library=variant.library_name, run_id=variant.run_id,
                                   chr=variant.chr, start=variant.pos, end=variant.end,
                                   gene=variant.gene, ref=variant.ref, alt=variant.alt, exon=variant.exon,
                                   codon=variant.codon_change, aa=variant.aa_change,
                                   rsids=",".join(variant.rs_ids),
                                   cosmic=",".join(variant.cosmic_ids) or None,
                                   cosmic_nsamples=variant.cosmic_data['num_samples'],
                                   cosmic_aa=variant.cosmic_data['aa'],
                                   amp=variant.amplicon_data['amplicon'],
                                   csig=variant.clinvar_data['significance'],
                                   cpath=variant.clinvar_data['pathogenic'],
                                   hgvs=variant.clinvar_data['hgvs'],
                                   cdis=variant.clinvar_data['disease'],
                                   biotype=variant.biotype, impact=variant.impact,
                                   impact_so=variant.impact_so, severity=variant.severity,
                                   max_maf_all=variant.max_maf_all,
                                   max_maf_no_fin=variant.max_maf_no_fin,
                                   max_som_aaf=variant.max_som_aaf,
                                   min_depth=variant.min_depth,
                                   max_depth=variant.max_depth,
                                   callers=",".join(variant.callers) or None))

            if 'mutect' in callers:
                report.write("\t{mfilter}\t{mmulti}\t{mdp}\t{mad}\t{maf}"
                             "".format(mfilter=variant.mutect.get('FILTER') or None,
                                       mmulti=variant.mutect.get('MULTIALLELIC') or None,
                                       mdp=variant.mutect.get('GTF_DP') or None,
                                       mad=variant.mutect.get('GTF_AD') or None,
                                       maf=variant.mutect.get('AAF') or None))

            if 'vardict' in callers:
                report.write("\t{vfilter}\t{vmulti}\t{vdp}\t{vad}\t{vaf}"
                             "".format(vfilter=variant.vardict.get('FILTER') or None,
                                       vmulti=variant.vardict.get('MULTIALLELIC') or None,
                                       vdp=variant.vardict.get('DP') or None,
                                       vad=variant.vardict.get('AD') or None,
                                       vaf=variant.vardict.get('AAF') or None))

            if 'freebayes' in callers:
                report.write("\t{ffilter}\t{fmulti}\t{fdp}\t{faf}\t{fro}\t{fao}"
                             "".format(ffilter=variant.freebayes.get('FILTER') or None,
                                       fmulti=variant.freebayes.get('MULTIALLELIC') or None,
                                       fdp=variant.freebayes.get('DP') or None,
                                       faf=variant.freebayes.get('AAF') or None,
                                       fro=variant.freebayes.get('RO') or None,
                                       fao=variant.freebayes.get('AO') or None))

            if 'scalpel' in callers:
                report.write("\t{sfilter}\t{smulti}\t{sdp}\t{sad}\t{saf}"
                             "".format(sfilter=variant.scalpel.get('FILTER') or None,
                                       smulti=variant.scalpel.get('MULTIALLELIC') or None,
                                       sdp=variant.scalpel.get('GTF_DP') or None,
                                       sad=variant.scalpel.get('GTF_AD') or None,
                                       saf=variant.scalpel.get('AAF') or None))

            if 'platypus' in callers:
                report.write("\t{plfilter}\t{plmulti}\t{pldp}\t{plad}\t{plaf}"
                             "".format(plfilter=variant.platypus.get('FILTER') or None,
                                       plmulti=variant.platypus.get('MULTIALLELIC') or None,
                                       plad=variant.platypus.get('TR') or None,
                                       pldp=variant.platypus.get('TC') or None,
                                       plaf=variant.platypus.get('AAF') or None))

            if 'pindel' in callers:
                report.write("\t{pfilter}\t{pmulti}\t{pdp}\t{pad}\t{paf}"
                             "".format(pfilter=variant.pindel.get('FILTER') or None,
                                       pmulti=variant.pindel.get('MULTIALLELIC') or None,
                                       pdp=variant.pindel.get('DP') or None,
                                       pad=variant.pindel.get('GTF_AD') or None,
                                       paf=variant.pindel.get('AAF') or None))

            report.write("\n")

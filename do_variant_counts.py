#!/usr/bin/env python
import sys
import getpass
import argparse
from ddb import configuration

from cassandra.cluster import Cluster
from cassandra.auth import PlainTextAuthProvider


def get_sample_variant_data(sample, samples, thresholds, authenticator):
    cluster = Cluster(['142.239.155.181', '142.239.155.182', '142.239.155.183',
                       '142.239.155.184'], auth_provider=authenticator)
    session = cluster.connect('variantstore')
    for library in samples[sample]:
        print samples[sample][library]['sample_name']
        rows = session.execute("""SELECT sample, run_id, reference_genome,
                               library_name, chr, pos, ref, alt, target_pool,
                               panel_name, sequencer, end, callers, type,
                               subtype, rs_id, rs_ids, cosmic_ids, gene,
                               transcript, exon, codon_change, aa_change,
                               biotype, severity, impact, impact_so, genes,
                               transcripts_data, in_cosmic, in_clinvar,
                               is_pathogenic, is_coding, is_lof, is_splicing,
                               population_freqs, clinvar_data, cosmic_data,
                               amplicon_data, max_maf_all, max_maf_no_fin,
                               min_depth, max_depth, min_som_aaf, max_som_aaf,
                               variant_filters, variant_categorization,
                               freebayes, mutect, scalpel, vardict, pindel,
                               platypus FROM sample_variant WHERE
                               sample=%s AND run_id=%s AND
                               reference_genome=%s AND library_name=%s""",
                               ([samples[sample][library]['sample_name'],
                                 samples[sample][library]['run_id'],
                                 config['genome_version'],
                                 samples[sample][library]['library_name']]))
        no_amplicon = 0
        num_rows = 0
        for variant_row in rows:
            num_rows += 1
            if variant_row.amplicon_data['amplicon'] == 'None':
                no_amplicon += 1
            else:
                print variant_row.sample,
                variant_row.amplicon_data['amplicon'],
                variant_row.chr,
                variant_row.pos,
                variant_row.ref,
                variant_row.alt
        print no_amplicon, num_rows
        print "Finished Variant Sample"
    print "Finished Variant Samples"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file',
                        help="Input configuration file for samples")
    parser.add_argument('-r', '--report',
                        help="Root name for reports (per sample)",
                        default='report')
    parser.add_argument('-a', '--address',
                        help="IP Address for Cassandra connection",
                        default='127.0.0.1')
    parser.add_argument('-u', '--username',
                        help='Cassandra username for login',
                        default=None)

    args = parser.parse_args()

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username,
                                              password=password)
    else:
        auth_provider = None

    callers = ("mutect", "platypus", "vardict", "scalpel", "freebayes",
               "pindel")

    sys.stdout.write("Processing samples\n")
    for sample in samples:
        get_sample_coverage_data(sample, samples, thresholds, auth_provider)
        get_sample_variant_data(sample, samples, thresholds, auth_provider)
    print "Finished Samples"

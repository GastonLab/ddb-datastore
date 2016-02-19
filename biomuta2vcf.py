import sys
import argparse

# This only works for BioMuta2. BioMuta3 appears to have a different format and order of columns, and is not documented
# with a header

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', help="BioMuta CSV input file")
    parser.add_argument('-o', '--outfile', help="Name of output VCF file")
    args = parser.parse_args()

    header = '##fileformat=VCFv4.1\n' \
             '##INFO=<ID=BM_Uniprot,Number=1,Type=String,Description="Uniprot ID">\n' \
             '##INFO=<ID=BM_Gene,Number=1,Type=String,Description="Gene name">\n' \
             '##INFO=<ID=BM_PolyPhen,Number=1,Type=String,Description="Polyphen score">\n' \
             '##INFO=<ID=BM_PMID,Number=1,Type=String,Description="Pubmed ID(s)">\n' \
             '##INFO=<ID=BM_cancertype,Number=1,Type=String,Description="Cancer Type">\n' \
             '##INFO=<ID=BM_source,Number=1,Type=String,Description="Source of information">\n' \
             '##INFO=<ID=BM_status,Number=1,Type=String,Description="Status of variant">\n' \
             '##INFO=<ID=BM_func,Number=1,Type=String,Description="Functional consequence of variant">\n' \
             '##CADDCOMMENT=<ID=comment,comment="{comment}">\n' \
             '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'

    with open(args.outfile, 'w') as outfile:
        outfile.write(header)
        with open(args.infile, 'r') as infile:
            infile.readline()
            for line in infile.readlines():
                info = line.split('\t')

                # There are some malformed lines that are not complete. Skip them
                if len(info) < 15:
                    # sys.stderr.write("WARNING: Error processing line: {}\n".format(line))
                    continue

                position = info[3]
                pos_sect = position.split(':')

                # Some lines don't have a position. Skip them
                if len(pos_sect) < 2:
                    # sys.stderr.write("WARNING: Error processing coordinates ({}) for line {}\n".format(position,
                    # line))
                    continue

                chrom = pos_sect[0]
                positions = pos_sect[1].split('-')

                outfile.write("{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t".format(chrom=chrom, pos=positions[0],
                                                                                  ref=info[5], alt=info[6]))

                outfile.write("BM_Uniprot={uni};BM_Gene={gene};BM_PolyPhen={poly};BM_PMID={pmid};"
                              "BM_cancertype={type};BM_source={source};BM_status={status};"
                              "BM_func={func}\n".format(uni=info[0], gene=info[1], poly=info[10], pmid=info[11],
                                                        type=info[12], source=info[13], status=info[14], func=info[15]))

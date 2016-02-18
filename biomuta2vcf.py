import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', help="BioMuta CSV input file")
    parser.add_argument('-o', '--outfile', help="Name of output VCF file")
    args = parser.parse_args()

    header = '##fileformat=VCFv4.1\n' \
             '##INFO=<ID=raw,Number=1,Type=Float,Description="raw cadd score">\n' \
             '##INFO=<ID=phred,Number=1,Type=Float,Description="phred-scaled cadd score">\n' \
             '##CADDCOMMENT=<ID=comment,comment="{comment}">\n' \
             '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'

    with open(args.outfile, 'w') as outfile:
        outfile.write(header)
        with open(args.infile, 'r') as infile:
            infile.readline()
            for line in infile.readlines():
                if len(info) < 15:
                    sys.stderr.write("WARNING: Error processing line: {}\n".format(line))
                info = line.split('\t')
                position = info[3]
                pos_sect = position.split(':')
                chrom = pos_sect[0]
                positions = pos_sect[1].split('-')
                outfile.write("{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t".format(chrom=chrom, pos=positions[0],
                                                                                  ref=info[5], alt=info[6]))

                outfile.write("BM_Uniprot={uni}\tBM_Gene={gene}\tBM_PolyPhen={poly}\tBM_PMID={pmid}\t"
                              "BM_cancertype={type}\tBM_source={source}\tBM_status={status}\t"
                              "BM_func={func}\n".format(uni=info[0], gene=info[1], poly=info[10], pmid=info[11],
                                                        type=info[12], source=info[13], status=info[14], func=info[15]))

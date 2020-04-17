import vcf
import sys
import json
import cyvcf2
import argparse

from cyvcf2 import VCF


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',
                        help="Input TSO500 VCF file")
    parser.add_argument('-a', '--annotation',
                        help="JSON formatted TSO500 annotations file")
    args = parser.parse_args()

    sys.stdout.write("Parsing TSO500 VCF\n")
    vcf = VCF(args.input)

    with open(args.annotation, 'r') as anno_file:
        anno_data = json.load(anno_file)

    for variant in vcf:
        # Do Something

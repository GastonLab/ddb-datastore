#!/usr/bin/env python

import sys
from cyvcf2 import VCF


if __name__ == "__main__":
    infile = sys.argv[1]
    # outfile = sys.argv[2]

    reader = VCF(infile)
    for record in reader:
        print record
        print record.gt_depths
        print record.gt_alt_depths

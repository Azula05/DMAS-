#!/usr/bin/env python3

import argparse
import os
import re



parser = argparse.ArgumentParser(description='give arguments to main snp script')
parser.add_argument("-i", nargs=1, required=True, help="user input file")
parser.add_argument('-u', nargs=1, required=True, help='SNP url', metavar='SNP url')

args = parser.parse_args()
SNP_url = args.u[0]

seq_ID = os.path.splitext(args.i[0])[0]

### Retrieve chrom, start and end info

input_file = open(args.i[0]).readline().rstrip()
seq, pos, SNP_pos = input_file.split('\t')

chrom, start, end = pos.rstrip().replace(":", "-").split("-")

# Retrieve SNPs from database

os.system("bigBedToBed " + SNP_url + " -chrom="+ chrom + " -start=" + start + " -end=" + end + " snps_" + seq_ID + ".bed") # TEMP SOLUTION, SEE BASECAMP

#open("snps_" + seq_ID + ".bed", 'w')


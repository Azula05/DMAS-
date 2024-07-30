#!/usr/bin/env python3
"""
This script is meant to get the SNPs from the database using the user input file.
Positions that are found in the common SNP database are avoided because they can introduce an extra mismatch.
"""
# Import the required libraries
import argparse
import os
import re

# Define the arguments that the user must provide
parser = argparse.ArgumentParser(description='give arguments to main snp script')
parser.add_argument("-i", nargs=1, required=True, help="user input file")
parser.add_argument('-u', nargs=1, required=True, help='SNP url', metavar='SNP url') 
## Example: http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.bb

args = parser.parse_args()
SNP_url = args.u[0]
############################################################################################################
#######################################   GET INFORMATION   ################################################
############################################################################################################
# Get the sequence ID from the input file
seq_ID = os.path.splitext(args.i[0])[0]
seq_ID = seq_ID.split("/")[-1]

# Retrieve chrom, start and end info
input_file = open(args.i[0]).readline().rstrip()
seq, pos, SNP_pos = input_file.split('\t')
chrom, start, end = pos.rstrip().replace(":", "-").split("-")

# Retrieve SNPs from database
os.system("bigBedToBed " + SNP_url + " -chrom="+ chrom + " -start=" + start + " -end=" + end + " snps_" + seq_ID + ".bed") # TEMP SOLUTION, SEE BASECAMP

#open("snps_" + seq_ID + ".bed", 'w')


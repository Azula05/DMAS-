#!/usr/bin/env python3
"""
This script will parse the Primer3 output file and select relevant information on the common primer.
It will then order the primers based on the melting temperature and amplicon length.
The theoretical best primer will be at the top based on these criteria.
"""
# ARGUMENT HANDLING
import argparse
import re
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import os

parser = argparse.ArgumentParser(description='give arguments to main snp script')
parser.add_argument("-i", nargs=1, required=True, help="Primer3 output file")
parser.add_argument("-w", nargs=1, required=True, help="Warning file")
parser.add_argument("-dnac", nargs=1, required=True, help="concetration of the oligos µM")
parser.add_argument("-Na", nargs=1, required=True, help="Na concentration mM")
parser.add_argument("-K", nargs=1, required=True, help="K concentration mM")
parser.add_argument("-Tris", nargs=1, required=True, help="Tris concentration mM")
parser.add_argument("-Mg", nargs=1, required=True, help="Mg concentration mM")
parser.add_argument("-dNTPs", nargs=1, required=True, help="dNTPs concentration mM")
args = parser.parse_args()

# files
common_primer_file = open(args.i[0], "r")
warning_file = open(args.w[0], "a")

# concentrations for temperature prediction
dnac = float(args.dnac[0])
Na_conc= float(args.Na[0])
K_conc = float(args.K[0])
Tris_conc = float(args.Tris[0])
Mg_conc = float(args.Mg[0])
dNTPs_conc = float(args.dNTPs[0])

####################################################################################################
##############################   Parse the Primer3 output file   ###################################
####################################################################################################
i = 0
for line in common_primer_file:
    # ID
    if "SEQUENCE_ID" in line:
        sequence_id = line.split("=")[1].strip().split("_")[0]
    # NR_of_PRIMERS_found
    if "PRIMER_PAIR_NUM_RETURNED" in line and "_SEQUENCE" in line:
        primers_returned = int(line.split("=")[1].strip())
        if primers_returned == 0:
            warning_file.write("No common primers found for sequence {}\n".format(sequence_id))
    pattern_L = r"PRIMER_LEFT_\d+_SEQUENCE="
    pattern_R = r"PRIMER_RIGHT_\d+_SEQUENCE="
    if re.search(pattern_L, line):
        left_primer = line.split("=")[1].strip()
        # MELTING TEMPERATURE LEFT
        myseq = Seq(left_primer)
        Tm_left = mt.Tm_NN(myseq, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
    if re.search(pattern_R, line):
        right_primer = line.split("=")[1].strip()
        # MELTING TEMPERATURE RIGHT
        myseq = Seq(right_primer)
        Tm_right = mt.Tm_NN(myseq, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
    pattern = r"PRIMER_LEFT_\d+_GC_PERCENT="
    if re.search(pattern, line):
        left_gc = line.split("=")[1].strip()
    pattern = r"PRIMER_RIGHT_\d+_GC_PERCENT="
    if re.search(pattern, line):
        right_gc = line.split("=")[1].strip()
    pattern = r"PRIMER_PAIR_\d+_PRODUCT_SIZE="
    if re.search(pattern, line):
        product_size = line.split("=")[1].strip()

    i += 1

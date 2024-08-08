#!/usr/bin/env python3
"""
This script will use user settings to create primer3 input files to find a common primer.
"""
# ARGUMENT HANDLING
import argparse
from Bio.Seq import Seq
import re
import os

parser = argparse.ArgumentParser(description='give arguments to main snp script')
parser.add_argument("-a", nargs=1, required=True, help="primer min left 3 prime distance") # Is dit nodig? wat doet dit?
parser.add_argument("-b", nargs=1, required=True, help="number of primers to return")
parser.add_argument("-c", nargs=1, required=True, help="min melting temperature")
parser.add_argument("-d", nargs=1, required=True, help="max melting temperature")
parser.add_argument("-e", nargs=1, required=True, help="optimal melting temperature")
parser.add_argument("-f", nargs=1, required=True, help="melting temperature difference between primers")
parser.add_argument("-g", nargs=1, required=True, help="min GC content")
parser.add_argument("-i", nargs=1, required=True, help="max GC content")
parser.add_argument("-j", nargs=1, required=True, help="optimal GC content")
parser.add_argument("-k", nargs=1, required=True, help="min amplicon length")
parser.add_argument("-l", nargs=1, required=True, help="max amplicon length")
parser.add_argument("-m", nargs=1, required=True, help="mispriming library")
parser.add_argument("-M", nargs=1, required=True, help="max mispriming library")

parser.add_argument("-t", nargs=1, required=True, help="seq-file with template")
parser.add_argument("-p", nargs=1, required=True, help="file with primers")

parser.add_argument("-dnac", nargs=1, required=True, help="concetration of the oligos µM")
parser.add_argument("-Na", nargs=1, required=True, help="Na concentration mM")
parser.add_argument("-K", nargs=1, required=True, help="K concentration mM")
parser.add_argument("-Tris", nargs=1, required=True, help="Tris concentration mM")
parser.add_argument("-Mg", nargs=1, required=True, help="Mg concentration mM")
parser.add_argument("-dNTPs", nargs=1, required=True, help="dNTPs concentration mM")

args = parser.parse_args()

# Set the user settings for the primer3 program.
MIN_LEFT_3_PRIME_DISTANCE = args.a[0]
NUMBER_OF_PRIMERS = args.b[0]
MIN_MELTING_TEMPERATURE = args.c[0]
MAX_MELTING_TEMPERATURE = args.d[0]
OPTIMAL_MELTING_TEMPERATURE = args.e[0]
MELTING_TEMPERATURE_DIFFERENCE = args.f[0]
MIN_GC_CONTENT = args.g[0]
MAX_GC_CONTENT = args.i[0]
OPTIMAL_GC_CONTENT = args.j[0]
MIN_AMPLICON_LENGTH = args.k[0]
MAX_AMPLICON_LENGTH = args.l[0]
MISPRIMING_LIBRARY = args.m[0]
MAX_MISPRIMING_LIBRARY = args.M[0]
input_file = args.t[0]
dnac_conc = args.dnac[0]
Na_conc = args.Na[0]
K_conc = args.K[0]
Tris_conc = args.Tris[0]
Mg_conc = args.Mg[0]
dNTPs_conc = args.dNTPs[0]

# get all the info from the sequence input file
seq_info = open(args.t[0]).readline().rstrip()
seq, pos, targeted_snp_pos = seq_info.split('\t')
chrom, start, end = pos.replace(':', '-').split('-')
forward = Seq(seq).split('[')[0]
wt = Seq(seq).split('[')[1].split('/')[0]
m = Seq(seq).split('[')[1].split('/')[1].split(']')[0]
reverse = (Seq(seq).split(']')[1])
template_seq = forward + wt + reverse
reverse = reverse.reverse_complement()

# get the sequence ID
seq_ID = os.path.splitext(args.t[0])[0]
seq_ID = seq_ID.split('/')[-1]

# get a forward and reverse primer from the primer file
primers = open(args.p[0], "r")
for line in primers:
    if re.search(r"DMAS-\d+_[A-Z0-9]+_F_WT", line):
        # primer3 needs an exact match so the primer is replaced by the exact match
        forward_len = len(line.split('\t')[1])
        template = forward + wt
        forward = template[-int(forward_len):]
        break
for line in primers:
    if re.search(r"DMAS-\d+_[A-Z0-9]+_R_WT", line):
        # primer3 needs an exact match so the primer is replaced by the exact match
        reverse_len = len(line.split('\t')[1])
        template = reverse + Seq(wt).complement()
        reverse = template[-int(reverse_len):]
        break
primers.close()

# write the forward input file for primer3
## settings in the input file will overwrite the settings in the primer3 settings file
Primer3_common = open("Primer3_" + seq_ID + "_common_primer_REV.txt", "w")
Primer3_common.write("SEQUENCE_ID=" + seq_ID + "_common_REV" "\n")
Primer3_common.write("SEQUENCE_TEMPLATE=" + str(template_seq) + "\n")
Primer3_common.write("SEQUENCE_PRIMER=" + str(forward) + "\n")
Primer3_common.write("PRIMER_NUM_RETURN=" + args.b[0] + "\n")
Primer3_common.write("PRIMER_MIN_THREE_PRIME_DISTANCE=" + args.a[0] + "\n")
Primer3_common.write("PRIMER_PRODUCT_SIZE_RANGE=" + args.k[0] + "-" + args.l[0] + "\n")
Primer3_common.write("PRIMER_MIN_TM=" + args.c[0] + "\n")
Primer3_common.write("PRIMER_MAX_TM=" + args.d[0] + "\n")
Primer3_common.write("PRIMER_OPT_TM=" + args.e[0] + "\n")
Primer3_common.write("PRIMER_PAIR_MAX_DIFF_TM=" + args.f[0] + "\n")
Primer3_common.write("PRIMER_MIN_GC=" + args.g[0] + "\n")
Primer3_common.write("PRIMER_MAX_GC=" + args.i[0] + "\n")
Primer3_common.write("PRIMER_OPT_GC_PERCENT=" + args.j[0] + "\n")
Primer3_common.write("PRIMER_OPT_GC_PERCENT=" + args.j[0] + "\n")
Primer3_common.write("PRIMER_MISPRIMING_LIBRARY=" + args.m[0] + "\n")
Primer3_common.write("PRIMER_INTERNAL_MISHYB_LIBRARY="+ args.m[0] + "\n")
Primer3_common.write("PRIMER_MAX_LIBRARY_MISPRIMING=" + args.M[0] + "\n")
Primer3_common.write("PRIMER_TASK=" + "generic" + "\n")
Primer3_common.write("PRIMER_PICK_LEFT_PRIMER=" + "0" + "\n")
Primer3_common.write("PRIMER_PICK_RIGHT_PRIMER=" + "1" + "\n")
Primer3_common.write("PRIMER_EXPLAIN_FLAG=" + "1" + "\n")
Primer3_common.write("PRIMER_DNA_CONC=" + dnac_conc +"\n")
Primer3_common.write("PRIMER_SALT_MONOVALENT="+ Na_conc + "\n")
Primer3_common.write("PRIMER_SALT_DIVALENT="+ Mg_conc + "\n")
Primer3_common.write("PRIMER_DNTP_CONC="+ dNTPs_conc + "\n")
Primer3_common.write("=\n")

Primer3_common.close()

# write the reverse input file for primer3
Primer3_common = open("Primer3_" + seq_ID + "_common_primer_FWD.txt", "w")
Primer3_common.write("SEQUENCE_ID=" + seq_ID + "_common_FWD" "\n")
Primer3_common.write("SEQUENCE_TEMPLATE=" + str(template_seq) + "\n")
Primer3_common.write("SEQUENCE_PRIMER_REVCOMP=" + str(reverse) + "\n")
Primer3_common.write("PRIMER_NUM_RETURN=" + args.b[0] + "\n")
Primer3_common.write("PRIMER_MIN_THREE_PRIME_DISTANCE=" + args.a[0] + "\n")
Primer3_common.write("PRIMER_PRODUCT_SIZE_RANGE=" + args.k[0] + "-" + args.l[0] + "\n")
Primer3_common.write("PRIMER_MIN_TM=" + args.c[0] + "\n")
Primer3_common.write("PRIMER_MAX_TM=" + args.d[0] + "\n")
Primer3_common.write("PRIMER_OPT_TM=" + args.e[0] + "\n")
Primer3_common.write("PRIMER_PAIR_MAX_DIFF_TM=" + args.f[0] + "\n")
Primer3_common.write("PRIMER_MIN_GC=" + args.g[0] + "\n")
Primer3_common.write("PRIMER_MAX_GC=" + args.i[0] + "\n")
Primer3_common.write("PRIMER_OPT_GC_PERCENT=" + args.j[0] + "\n")
Primer3_common.write("PRIMER_OPT_GC_PERCENT=" + args.j[0] + "\n")
Primer3_common.write("PRIMER_MISPRIMING_LIBRARY=" + args.m[0] + "\n")
Primer3_common.write("PRIMER_INTERNAL_MISHYB_LIBRARY="+ args.m[0] + "\n")
Primer3_common.write("PRIMER_MAX_LIBRARY_MISPRIMING=" + args.M[0] + "\n")
Primer3_common.write("PRIMER_TASK=" + "generic" + "\n")
Primer3_common.write("PRIMER_PICK_LEFT_PRIMER=" + "1" + "\n")
Primer3_common.write("PRIMER_PICK_RIGHT_PRIMER=" + "0" + "\n")
Primer3_common.write("PRIMER_EXPLAIN_FLAG=" + "1" + "\n")
Primer3_common.write("PRIMER_DNA_CONC=" + dnac_conc +"\n")
Primer3_common.write("PRIMER_SALT_MONOVALENT="+ Na_conc + "\n")
Primer3_common.write("PRIMER_SALT_DIVALENT="+ Mg_conc + "\n")
Primer3_common.write("PRIMER_DNTP_CONC="+ dNTPs_conc + "\n")
Primer3_common.write("=\n")

Primer3_common.close()
#!/usr/bin/env python3
"""
This script will use user settings to create an upfront filter.
This upfront filter will use information from the secondary structure and SNP files to avoid certain regions.
Based on this information a Primer3 input file will be created to design common primers.
"""
# ARGUMENT HANDLING
import argparse
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
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
parser.add_argument("-q", nargs=1, required=True, help="variable for filtering options: no/yes/snp/str")
parser.add_argument("-s", nargs=1, required=True, help="bed file with SNPs")
parser.add_argument("-t", nargs=1, required=True, help="seq-file with template")
parser.add_argument("-u", nargs=1, required=True, help="file with secondary structures")
parser.add_argument("-p", nargs=1, required=True, help="Position of the mismatch eighter all, 2, 3 or 4 postions to 5' end")
parser.add_argument("-dnac", nargs=1, required=True, help="concetration of the oligos µM")
parser.add_argument("-Na", nargs=1, required=True, help="Na concentration mM")
parser.add_argument("-K", nargs=1, required=True, help="K concentration mM")
parser.add_argument("-Tris", nargs=1, required=True, help="Tris concentration mM")
parser.add_argument("-Mg", nargs=1, required=True, help="Mg concentration mM")
parser.add_argument("-dNTPs", nargs=1, required=True, help="dNTPs concentration mM")

args = parser.parse_args()

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
filter_on = args.q[0]
SNP_file = args.s[0]
seq_file = args.t[0]
sec_str_file = args.u[0]
position_mismatch = args.p[0]
dnac = args.dnac[0]
Na = args.Na[0]
K = args.K[0]
Tris = args.Tris[0]
Mg = args.Mg[0]
dNTPs = args.dNTPs[0]

# Create a input file for the primer3_core program to design common primers.
## Uses the whole template sequence to design primers.
## The position of the SNP needs to be excluded otherwise it is not common.

# get all the info from the sequence input file
seq_info = open(args.t[0]).readline().rstrip()
seq, pos, targeted_snp_pos = seq_info.split('\t')
chrom, start, end = pos.replace(':', '-').split('-')
start = int(start)
end = int(end)
targeted_snp_pos = int(targeted_snp_pos)+1 # 1-based
template_seq = Seq(seq).split('[')[0] + Seq(seq).split('[')[1].split('/')[0] + Seq(seq).split(']')[1]
seq_length = len(template_seq) # 1-based
seq_ID = os.path.splitext(args.t[0])[0]
seq_ID = seq_ID.split('/')[-1]

Primer3_common = open("Primer3_common_primers_"+ seq_ID +".txt", "w")
Primer3_common.write("SEQUENCE_ID=" + seq_ID + "_common" "\n")
Primer3_common.write("SEQUENCE_TEMPLATE=" + str(template_seq) + "\n")
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

####################################################################################################
#####################################   structure filter   #########################################
####################################################################################################
avoid_range = []
# if this filter is on (yes meaning both str and snp)
if filter_on == 'yes' or filter_on == "str":
	folding_file = open(sec_str_file)

	# get folding info
	for line in folding_file:
		# check if line is not empty
		if line.rstrip() != '':
			ID, structure, deltaG, fold_temp_avoid = line.rstrip().split('\t')
		# put folding info into format for primer 3 (space-separated: position,length)
		## only if there are positions to avoid
		if fold_temp_avoid != '[]':
			fold_temp_avoid = fold_temp_avoid.replace('[', "").replace(']', '').split(', ')
			fold_temp_avoid = [int(x) for x in fold_temp_avoid]

			begin_region = fold_temp_avoid[0]
			index = 0
			length = 1
			# go over all positions and add them to avoid_range
			for position in fold_temp_avoid:
				if index + 1 < len(fold_temp_avoid):
					if position == fold_temp_avoid[index + 1] - 1:
						length += 1
						index += 1
					else:
						avoid_range.append([begin_region, length])
						begin_region = fold_temp_avoid[index + 1]
						index += 1
						length = 1
				else:
					avoid_range.append([begin_region, length])
					length = 1
					
####################################################################################################
####################################     SNP filter    #############################################
####################################################################################################
if filter_on == 'yes' or filter_on == "snp":
	# get SNP info
	SNPs = open(SNP_file)
	for line in SNPs:
		# relative position
		relative_position = int(line.split('\t')[1]) - (start) +2 # Make it 1-based (start position is 0-based +1, substraction needs +1)
		# length of the SNP
		length_SNP = int(line.split('\t')[2]) - int(line.split('\t')[1])
		avoid_range.append([relative_position, length_SNP])
	# get all SNPs in the region
	
####################################################################################################
####################################    Write positions to avoid    ################################
####################################################################################################
if filter_on != 'no':
    avoid_range.append([targeted_snp_pos, "1"]) # make sure the SNP is avoided, you do not want this in the commmon pirmer
    avoid_range = sorted(avoid_range, key=lambda x: x[0])
    avoid_str = ""
    for item in avoid_range:
        avoid_str = avoid_str + str(item[0]) + "," + str(item[1]) + " "
    avoid_str = avoid_str[:-1]
    Primer3_common.write("SEQUENCE_EXCLUDED_REGION=" + avoid_str + "\n")


# add last line with '=' again
Primer3_common.write("=\n")

# close all files
Primer3_common.close()
Primer3_common.close()
#!/usr/bin/env python3
"""
This script will add pass or fail to the table based on the filters used.
These tags can be used to filter the tabel once more guideliness have been set.
The tags are used to create a log file as well which reports how many primers passed or failed each filter.
The output of this will include the full table a log file and a filtered table.
"""

import argparse
import os

parser = argparse.ArgumentParser(description="give arguments to filter script")
parser.add_argument("-i", nargs = 1, required=True, help="file with sequence info example: DMAS-0.txt")
parser.add_argument("-p", nargs = 1, required=True, help="file with all designed primers")
parser.add_argument("-s", nargs = 1, required=True, help="specificity filter: off, strict, loose")
parser.add_argument("-S", nargs=1, required=True, help="SNP filter: off, strict, loose")
parser.add_argument("-t", nargs=1, required=True, help="Secondary structure filter: off, strict, loose")
parser.add_argument("-v", nargs=1, required=True, help="Validation filter off, strict, loose") 
args = parser.parse_args()

# get the arguments
input_file = args.i[0]
primers_file = args.p[0]
specificity_filter = args.s[0]
SNP_filter = args.S[0]
secondary_structure_filter = args.t[0]
validation_filter = args.v[0]

# Open the table and save the contents
with open(primers_file,'r') as table:
    lines = table.readlines()
table.close()

# write a new table
full_table = open(primers_file, 'w')
# get the sequence ID
seq_ID = os.path.splitext(args.i[0])[0]
seq_ID = seq_ID.split("/")[-1]
# create a filtered table
filtered_table = open(seq_ID + "_filtered.tsv", 'w')

# go through the table
line_nr = 0
for line in lines:
    line = line.strip().split('\t')
    # Ignore the header and write a new one
    if line_nr == 0:
        full_table.write("Name\tSpecific_primer\tMatch_Tm\tSingle_MM_Tm\tDouble_MM_Tm\tMM_delta\tGC%\tLenght\tCommon_primer\tMatch_Tm_common\tGC%_common\tLength_common\tSNPs_FWD\tSec_str_FWD\tSNPs_REV\tSec_str_REV\tDeltaG_template\tFWD_validation\tREV_validation\tAmplicon\tAmp_length\tPredicted_structure\tAmplicon_delta_G\tForward_specificity\tReverse_specificity\tSpecificity_filter\tSNP_filter\tSec_str_filter\tValidation_filter\n")
        line_nr += 1
        continue
    ################################################################################################
    ###################################   Specificity   ############################################
    ################################################################################################
    # rules from https://academic.oup.com/clinchem/article/59/10/1470/5622018?login=true fig 6 are used as filter criteria
    # The criteria are:
    # Strict: 
    # - no mismatches in either primer = off-target
    # - primers with at least 4 mismatches for a single primer = no off-target
    # - primers with a total of at least 5 mismatches between both primers = no off-target
    # Loose:
    # - primers with at least 3 mismatches for a single primer = no off-target
    # - primers with a total of at least 4 mismatches between both primers = no off-target

    # MM primer1	MM primer2	sum MM	loose			strict
    # 0				0			0		off-target		off-target
    # 1				0			1		off-target		off-target
    # 0				1			1		off-target		off-target
    # 2				0			2		off-target		off-target
    # 0				2			2		off-target		off-target
    # 1				1			2		off-target		off-target
    # 2				1			3		off-target		off-target
    # 1				2			3		off-target		off-target
    # 3				0			3		no off-target	off-target
    # 0				3			3		no off-target	off-target
    # 1				3			4		no off-target	off-target
    # 3				1			4		no off-target	off-target
    # 2				2			4		no off-target	off-target
    # 4				0			4		no off-target	no off-target
    # 0				4			4		no off-target	no off-target
    # 2				3			5		no off-target	no off-target
    # 3				2			5		no off-target	no off-target
    # 4				1			5		no off-target	no off-target
    # 1				4			5		no off-target	no off-target
    # 3				3			6		no off-target	no off-target

    specifity_lost = 0
    if specificity_filter == "off":
        Specificity_tag = "NA"
    elif specificity_filter == "loose":
        try:
            Forward_specificity = line[23]
            i=0
            Forward_specificity = Forward_specificity.replace("[","").replace("]","").replace("'","")
            Forward_specificity = Forward_specificity.split(":")
            print(Forward_specificity)
            Reverse_specificity = line[24]
        except:
            specificity_tag = "error"



    ################################################################################################
    ###################################   SNPs   ##################################################
    ################################################################################################

    ################################################################################################
    ###################################   Secondary structure   ####################################
    ################################################################################################

    ################################################################################################
    ###################################   Validation   #############################################
    ################################################################################################
    line_nr += 1

total = line_nr # 1-based counting (0 is header)
full_table.close()
filtered_table.close()

################################################################################################
# Append to the log file
################################################################################################
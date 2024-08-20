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
    Specificity_tag = "PASS"
    # filter off => write NA
    if specificity_filter != "off":
        try:
            # get the number of mismatches
            Forward_specificity = line[23]
            Forward_specificity = Forward_specificity.replace("[","").replace("]","").replace("'","")
            Forward_specificity = Forward_specificity.split(":")
            Reverse_specificity = line[24]
            Reverse_specificity = Reverse_specificity.replace("[","").replace("]","").replace("'","")
            Reverse_specificity = Reverse_specificity.split(":")
            # get the original input coordinates
            input_file = open(input_file, 'r').readline()
            coord = input_file.split("\t")[-2]
            input_file.close()
        except:
            specificity_tag = "error"

    # filter loose
    def specificity_loose(nm_forward, nm_reverse):
        if nm_forward < 3 and nm_reverse < 3:
            if nm_forward == 2 and nm_reverse == 2:
                Specificity_tag = "PASS"
            else:
                Specificity_tag = "FAIL"
        else:
            Specificity_tag = "PASS"
        return Specificity_tag



    if specificity_filter == "loose":
        # placeholder values
        nm_forward_prev = 10
        nm_reverse_prev = 10
        # The find the lowest amount of mismatches in off-targets: if to low fail
        for i in range(0,10):
            try:
                chr_forward = Forward_specificity[i]
                start_forward = Forward_specificity[i+1]
                nm_forward = int(Forward_specificity[i+3])
            except:
                continue

            try:
                chr_reverse = Reverse_specificity[i]
                start_reverse = Reverse_specificity[i+1]
                nm_reverse = int(Reverse_specificity[i+3])
            except:
                continue
            # Once a primer pair has failed stop checking the rest
            if Specificity_tag == "FAIL":
                continue
            # check if not on target:
            if chr_forward == coord.split(":")[0] and chr_reverse == coord.split(":")[0]:
                if int(start_forward) > int(coord.split(":")[1].split("-")[0]) and int(start_forward) < int(coord.split(":")[1].split("-")[1]) and int(start_reverse) > int(coord.split(":")[1].split("-")[0]) and int(start_reverse) < int(coord.split(":")[1].split("-")[1]):
                    Specificity_tag = "PASS"
                # On the same chromosome but not the same region
                else:
                    if nm_forward_prev < nm_forward:
                        nm_forward = nm_forward_prev
                    elif nm_reverse_prev < nm_reverse:
                        nm_reverse = nm_reverse_prev
                    Specificity_tag = specificity_loose(nm_forward, nm_reverse)
            # off target
            else:
                if nm_forward_prev < nm_forward:
                    nm_forward = nm_forward_prev
                elif nm_reverse_prev < nm_reverse:
                    nm_reverse = nm_reverse_prev
                Specificity_tag = specificity_loose(nm_forward, nm_reverse)


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

    ################################################################################################
    ###################################   Write the table   ########################################
    ################################################################################################
    full_table.write("\t".join(line) + "\t" + Specificity_tag + "\n")

full_table.close()
filtered_table.close()

################################################################################################
# Append to the log file
################################################################################################
total = line_nr # 1-based counting (0 is header)
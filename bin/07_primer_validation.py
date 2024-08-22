#!/usr/bin/env python3
"""
This script will add information about the primers to the csv file with the designed primers.
- It will add the best common primer to the table as well as write the alternatives to a separate file.
- It will use primer3_core to validate the primers.
- It will look if any snips were found in the primers.
- It will check if any secondary structures were found in the primers.
- It will add the amplicon to the table.
"""

import argparse
import re
import os
from Bio.Seq import Seq
import pandas as pd

parser = argparse.ArgumentParser(description='give arguments to main primer validation script')
parser.add_argument("-F", nargs=1, required=True, help="file with the forward common primers output")
parser.add_argument("-R", nargs=1, required=True, help="file with the reverse common primers output")
parser.add_argument("-T", nargs=1, required=True, help="file with the template sequence; input file")
parser.add_argument("-s", nargs=1, required=True, help="file with the SNP information")
parser.add_argument("-u", nargs=1, required=True, help="file with the secondary structure information for the wild type")
parser.add_argument("-U", nargs=1, required=True, help="file with the secondary structure information for the mutant")
parser.add_argument("-o", nargs=1, required=True, help="tsv file with the primer information")
parser.add_argument("-p", nargs=1, required=True, help="file with the primer3 settings")
parser.add_argument("-q", nargs=1, required=True, help="upfront checks: yes, snp, str, no")
args = parser.parse_args()

# Get the file names from the arguments
file_forward = args.F[0]
file_reverse = args.R[0]
file_template = args.T[0]
file_snps = args.s[0]
file_structures_wt = args.u[0]
file_structures_mut = args.U[0]
file_output = args.o[0]
primer3_settings= args.p[0]
# get the ID from the file
ID = file_template.split("_")[0].replace(".txt", "")
ID = ID.split("/")[-1]
# upfront check
checks = args.q[0]
SNP_avoid_range = {}
####################################################################################################
##############################   Add common primers to the table   #################################
####################################################################################################

# Start a file to write all alternatives for the common primers to
alternatives = open(ID + "_common_alternatives.txt", 'w')
alternatives.write("Common FORWARD primer alternatives\n")

# check all the common forward primers
for i in range(0, 21): # 1-20
    Forward = open(file_forward, 'r')
    for line in Forward:
        # get the sequence template
        pattern = r"SEQUENCE_TEMPLATE"
        if re.search(pattern, line):
            template = line.split("=")[1].strip()
        # get the primer sequence
        pattern = r"PRIMER_LEFT_" + str(i) + "_SEQUENCE="
        if re.search(pattern, line):
            primer = line.split("=")[1].strip()
            num = int(line.split("_")[2])
        # get the primer melting temperature
        pattern = r"PRIMER_LEFT_" + str(i) + "_TM="
        if re.search(pattern, line):
            tm = line.split("=")[1].strip()
        # get the primer GC content
        pattern = r"PRIMER_LEFT_" + str(i) + "_GC_PERCENT="
        if re.search(pattern, line):
            gc = line.split("=")[1].strip()
    # If there are less than 20 primers, stop the loop earlier
    if num != i: 
        break
    # The first primer is the best common primer and will be used in the table
    elif num == 0:
        common_FWD = primer
        common_TM = tm
        common_GC = gc
    # back-up primer
    elif num == 1:
        common_FWD_1 = primer
        common_TM_1 = tm
        common_GC_1 = gc
        alternatives.write(ID + "_common_F_" + str(num) + "\t" + primer + "\t" + tm + "\t" + gc + "\n")
    # The other primers are alternatives and will be written to the alternatives
    elif num > 1:
        alternatives.write(ID + "_common_F_" + str(num) + "\t" + primer + "\t" + tm + "\t" + gc + "\n")
Forward.close()

# check all common reverse primers
alternatives.write("\nCommon REVERSE primer alternatives\n")

for i in range(0, 21):
    Reverse = open(file_reverse, 'r')
    for line in Reverse: 
        # get the primer sequence
        pattern = r"PRIMER_RIGHT_" + str(i) + "_SEQUENCE="
        if re.search(pattern, line):
            primer = line.split("=")[1].strip()
            num = int(line.split("_")[2])
        # get the primer melting temperature
        pattern = r"PRIMER_RIGHT_" + str(i) + "_TM="
        if re.search(pattern, line):
            tm = line.split("=")[1].strip()
        # get the primer GC content
        pattern = r"PRIMER_RIGHT_" + str(i) + "_GC_PERCENT="
        if re.search(pattern, line):
            gc = line.split("=")[1].strip()
    # If there are less than 20 primers, stop the loop earlier
    if num != i:
        break
    # The first primer is the best common primer and will be used in the table
    elif num == 0:
        common_REV = primer
        common_TM_REV = tm
        common_GC_REV = gc
    elif num == 1:
        common_REV_1 = primer
        common_TM_REV_1 = tm
        common_GC_REV_1 = gc
        alternatives.write(ID + "_common_R_" + str(num) + "\t" + primer + "\t" + tm + "\t" + gc + "\n")
    # The other primers are alternatives and will be written to the alternatives
    elif num > 1:
        alternatives.write(ID + "_common_R_" + str(num) + "\t" + primer + "\t" + tm + "\t" + gc + "\n")
# close the files
Forward.close()
alternatives.close()
Reverse.close()

####################################################################################################
######################################   Validate primers   ########################################
####################################################################################################

def template_adaptation(left, right, name):
    change = name.split("_")[1]
    if change == "2A":
        num = 2
        change = "A"
    elif change == "2G":
        num = 2
        change = "G"
    elif change == "2T":
        num = 2
        change = "T"
    elif change == "2C":
        num = 2
        change = "C"
    elif change == "3A":
        num = 3
        change = "A"
    elif change == "3G":
        num = 3
        change = "G"
    elif change == "3T":
        num = 3
        change = "T"
    elif change == "3C":
        num = 3
        change = "C"
    elif change == "4A":
        num = 4
        change = "A"
    elif change == "4G":
        num = 4
        change = "G"
    elif change == "4T":
        num = 4
        change = "T"
    elif change == "4C":
        num = 4
        change = "C"
    left = list(left)
    left[-num] = change
    left = "".join(left)
    right = str(Seq(right).reverse_complement())
    right = list(right)
    right[-num] = change
    right = "".join(right)
    right = str(Seq(right).reverse_complement())
    return left , right

# save the old lines in a list
with open(file_output, 'r') as table:
    lines = table.readlines()
template_file = open(file_template, 'r')
# Info from the input file
for line in template_file:
    if line != "":
        seq, pos, position_of_interest = line.split('\t')
        chrom, start, end = pos.replace(':', '-').split('-')
        start = int(start)
        end = int(end)
        position_of_interest = int(position_of_interest)
        seq_length = len(seq)-4
        left = line.split("[")[0]
        wt = line.split("[")[1].split("]")[0].split("/")[0]
        m = line.split("]")[0].split("/")[1]
        right = line.split("]")[1].split("\t")[0]

# append the common priemers to the line and write the new line to the file
with open(file_output,'w') as output:
    # Header line
    output.write("Name\tSpecific_primer\tMatch_Tm\tSingle_MM_Tm\tDouble_MM_Tm\tMM_delta\tGC%\tLenght\tCommon_primer\tMatch_Tm_common\tGC%_common\tLength\tFWD_validation\tREV_validation\n")
    for line in lines:
        # Forward primers
        if "_F_WT" in line:
            # validate the primers
            line_list = line.split("\t")
            name = line_list[0]
            changed_L, changed_R = template_adaptation(left, right ,name)
            template = changed_L + wt + right
            forward = line_list[1]
            # validation
            command = "SEQUENCE_ID=" + name + "\nSEQUENCE_TEMPLATE=" + template + "\nSEQUENCE_PRIMER=" + forward + "\nSEQUENCE_PRIMER_REVCOMP=" + common_REV + "\nPRIMER_TASK=check_primers" + "\nPRIMER_EXPLAIN_FLAG=1 \nPRIMER_MIN_TM=45 \nPRIMER_MAX_TM=63 \nPRIMER_MIN_SIZE=16 \nPRIMER_MAX_SIZE=36 \n="
            try:
                process = os.popen("echo \"" + command + "\" | primer3_core -p3_settings_file=" + primer3_settings)            
                validation = process.read()
                process.close()  # Ensure proper resource management
                # Extract stdout from the command
                arguments = validation.split("\n")
                left_validation = arguments[10].split("=")[1]
                right_validation = arguments[11].split("=")[1]
                # change the line
                line = line.rstrip() + "\t" + common_REV + "\t" + common_TM_REV + "\t" + common_GC_REV + "\t" + str(len(common_REV)) + "\t" + left_validation + "\t" + right_validation + "\n"
                output.write(line)
            except:
                # try backup primer
                try:
                    back_up_F_WT = "yes"
                    command = "SEQUENCE_ID=" + name + "\nSEQUENCE_TEMPLATE=" + template + "\nSEQUENCE_PRIMER=" + forward + "\nSEQUENCE_PRIMER_REVCOMP=" + common_REV_1 + "\nPRIMER_TASK=check_primers" + "\nPRIMER_EXPLAIN_FLAG=1 \nPRIMER_MIN_TM=45 \nPRIMER_MAX_TM=63 \nPRIMER_MIN_SIZE=16 \nPRIMER_MAX_SIZE=36 \n="
                    process = os.popen("echo \"" + command + "\" | primer3_core -p3_settings_file=" + primer3_settings)            
                    validation = process.read()
                    process.close()  # Ensure proper resource management
                    # Extract stdout from the command
                    arguments = validation.split("\n")
                    left_validation = arguments[10].split("=")[1]
                    right_validation = arguments[11].split("=")[1]
                    # change the line
                    line = line.rstrip() + "\t" + common_REV_1 + "\t" + common_TM_REV_1 + "\t" + common_GC_REV_1 + "\t" + str(len(common_REV_1)) + "\t" + left_validation + "\t" + right_validation + "\n"
                    output.write(line)
                except:
                    line = line.rstrip() + "\t" + common_REV + "\t" + common_TM_REV + "\t" + common_GC_REV + "\t" + str(len(common_REV)) + "\t" + "primer3 failed to validate primers" + "\t" + "primer3 failed to validate primers"+ "\n"
                    output.write(line)

        if "_F_MUT" in line:
            # validate the primers
            line_list = line.split("\t")
            name = line_list[0]
            ## needs to be an exact match
            changed_L, changed_R = template_adaptation(left, right ,name)
            template_MUT = changed_L + m + right
            forward = line_list[1]
            # validation
            command = "SEQUENCE_ID=" + name + "\nSEQUENCE_TEMPLATE=" + template_MUT + "\nSEQUENCE_PRIMER=" + forward + "\nSEQUENCE_PRIMER_REVCOMP=" + common_REV + "\nPRIMER_TASK=check_primers" + "\nPRIMER_EXPLAIN_FLAG=1 \nPRIMER_MIN_TM=45 \nPRIMER_MAX_TM=63 \nPRIMER_MIN_SIZE=16 \nPRIMER_MAX_SIZE=36 \n="
            try:
                process = os.popen("echo \"" + command + "\" | primer3_core -p3_settings_file=" + primer3_settings)            
                validation = process.read()
                process.close()  # Ensure proper resource management
                # Extract stdout from the command
                arguments = validation.split("\n")
                left_validation = arguments[10].split("=")[1]
                right_validation = arguments[11].split("=")[1]
                # change the line
                line = line.rstrip() + "\t" + common_REV + "\t" + common_TM_REV + "\t" + common_GC_REV + "\t" + str(len(common_REV)) + "\t" + left_validation + "\t" + right_validation + "\n"
                output.write(line)
            except:
                try:
                    back_up_F_MUT = "yes"
                    command = "SEQUENCE_ID=" + name + "\nSEQUENCE_TEMPLATE=" + template_MUT + "\nSEQUENCE_PRIMER=" + forward + "\nSEQUENCE_PRIMER_REVCOMP=" + common_REV_1 + "\nPRIMER_TASK=check_primers" + "\nPRIMER_EXPLAIN_FLAG=1 \nPRIMER_MIN_TM=45 \nPRIMER_MAX_TM=63 \nPRIMER_MIN_SIZE=16 \nPRIMER_MAX_SIZE=36 \n="
                    process = os.popen("echo \"" + command + "\" | primer3_core -p3_settings_file=" + primer3_settings)            
                    validation = process.read()
                    process.close()  # Ensure proper resource management
                    # Extract stdout from the command
                    arguments = validation.split("\n")
                    left_validation = arguments[10].split("=")[1]
                    right_validation = arguments[11].split("=")[1]
                    # change the line
                    line = line.rstrip() + "\t" + common_REV_1 + "\t" + common_TM_REV_1 + "\t" + common_GC_REV_1 + "\t" + str(len(common_REV_1)) + "\t" + left_validation + "\t" + right_validation + "\n"
                    output.write(line)
                except:
                    line = line.rstrip() + "\t" + common_REV + "\t" + common_TM_REV + "\t" + common_GC_REV + "\t" + str(len(common_REV)) + "\t" + "primer3 failed to validate primers" + "\t" + "primer3 failed to validate primers"+ "\n"
                    output.write(line)
  
        # Reverse primers
        if "_R_WT" in line:
            # validate the primers
            line_list = line.split("\t")
            name = line_list[0]
            ## needs to be an exact match
            changed_L, changed_R = template_adaptation(left, right ,name)
            template = left + wt + changed_R
            reverse = line_list[1]
            # validation
            command = "SEQUENCE_ID=" + name + "\nSEQUENCE_TEMPLATE=" + template + "\nSEQUENCE_PRIMER=" + common_FWD + "\nSEQUENCE_PRIMER_REVCOMP=" + reverse + "\nPRIMER_TASK=check_primers" + "\nPRIMER_EXPLAIN_FLAG=1 \nPRIMER_MIN_TM=45 \nPRIMER_MAX_TM=63 \nPRIMER_MIN_SIZE=16 \nPRIMER_MAX_SIZE=36 \n="
            try:
                process = os.popen("echo \"" + command + "\" | primer3_core -p3_settings_file=" + primer3_settings)            
                validation = process.read()
                process.close()  # Ensure proper resource management
                # Extract stdout from the command
                arguments = validation.split("\n")
                left_validation = arguments[10].split("=")[1]
                right_validation = arguments[11].split("=")[1]
                # change the line
                line = line.rstrip() + "\t" + common_FWD + "\t" + common_TM + "\t" + common_GC + "\t" + str(len(common_FWD)) + "\t" + left_validation + "\t" + right_validation + "\n"
                output.write(line)
            except:
                try:
                    back_up_R_WT = "yes"
                    command = "SEQUENCE_ID=" + name + "\nSEQUENCE_TEMPLATE=" + template + "\nSEQUENCE_PRIMER=" + common_FWD_1 + "\nSEQUENCE_PRIMER_REVCOMP=" + reverse + "\nPRIMER_TASK=check_primers" + "\nPRIMER_EXPLAIN_FLAG=1 \nPRIMER_MIN_TM=45 \nPRIMER_MAX_TM=63 \nPRIMER_MIN_SIZE=16 \nPRIMER_MAX_SIZE=36 \n="
                    process = os.popen("echo \"" + command + "\" | primer3_core -p3_settings_file=" + primer3_settings)            
                    validation = process.read()
                    process.close()  # Ensure proper resource management
                    # Extract stdout from the command
                    arguments = validation.split("\n")
                    left_validation = arguments[10].split("=")[1]
                    right_validation = arguments[11].split("=")[1]
                    # change the line
                    line = line.rstrip() + "\t" + common_FWD_1 + "\t" + common_TM_1 + "\t" + common_GC_1 + "\t" + str(len(common_FWD_1)) + "\t" + left_validation + "\t" + right_validation + "\n"
                    output.write(line)
                except:
                    line = line.rstrip() + "\t" + common_FWD + "\t" + common_TM + "\t" + common_GC + "\t" + str(len(common_FWD)) + "\t" + "primer3 failed to validate primers" + "\t" + "primer3 failed to validate primers"+ "\n"
                    output.write(line)
        if "_R_MUT" in line:
            # validate the primers
            line_list = line.split("\t")
            name = line_list[0]
            ## needs to be an exact match
            changed_L, changed_R = template_adaptation(left, right ,name)
            template_MUT = left + m + changed_R
            reverse = line_list[1]
            # validation
            command = "SEQUENCE_ID=" + name + "\nSEQUENCE_TEMPLATE=" + template_MUT + "\nSEQUENCE_PRIMER=" + common_FWD + "\nSEQUENCE_PRIMER_REVCOMP=" + reverse + "\nPRIMER_TASK=check_primers" + "\nPRIMER_EXPLAIN_FLAG=1 \nPRIMER_MIN_TM=45 \nPRIMER_MAX_TM=63 \nPRIMER_MIN_SIZE=16 \nPRIMER_MAX_SIZE=36 \n="
            try:
                process = os.popen("echo \"" + command + "\" | primer3_core -p3_settings_file=" + primer3_settings)            
                validation = process.read()
                process.close()  # Ensure proper resource management
                # Extract stdout from the command
                arguments = validation.split("\n")
                left_validation = arguments[10].split("=")[1]
                right_validation = arguments[11].split("=")[1]
                # change the line
                line = line.rstrip() + "\t" + common_FWD + "\t" + common_TM + "\t" + common_GC + "\t" + str(len(common_FWD)) + "\t" + left_validation + "\t" + right_validation + "\n"
                output.write(line)
            except:
                try:
                    back_up_R_MUT = "yes"
                    command = "SEQUENCE_ID=" + name + "\nSEQUENCE_TEMPLATE=" + template_MUT + "\nSEQUENCE_PRIMER=" + common_FWD_1 + "\nSEQUENCE_PRIMER_REVCOMP=" + reverse + "\nPRIMER_TASK=check_primers" + "\nPRIMER_EXPLAIN_FLAG=1 \nPRIMER_MIN_TM=45 \nPRIMER_MAX_TM=63 \nPRIMER_MIN_SIZE=16 \nPRIMER_MAX_SIZE=36 \n="
                    process = os.popen("echo \"" + command + "\" | primer3_core -p3_settings_file=" + primer3_settings)            
                    validation = process.read()
                    process.close()  # Ensure proper resource management
                    # Extract stdout from the command
                    arguments = validation.split("\n")
                    left_validation = arguments[10].split("=")[1]
                    right_validation = arguments[11].split("=")[1]
                    # change the line
                    line = line.rstrip() + "\t" + common_FWD_1 + "\t" + common_TM_1 + "\t" + common_GC_1 + "\t" + str(len(common_FWD_1)) + "\t" + left_validation + "\t" + right_validation + "\n"
                    output.write(line)
                except:
                    line = line.rstrip() + "\t" + common_FWD + "\t" + common_TM + "\t" + common_GC + "\t" + str(len(common_FWD))  + "\t" + "primer3 failed to validate primers" + "\t" + "primer3 failed to validate primers"+ "\n"
                    output.write(line)
output.close()
table.close()

####################################################################################################
#######################################   Check for SNPs   #########################################
####################################################################################################
if checks == "yes" or checks =="YES" or checks == "snp" or checks == "SNP":
    # open the SNP file
    SNPs = open(args.s[0])
    # loop all the SNPs that where found
    for line in SNPs:
        # if the line is not empty
        if line != "":
            line = line.split("\t")
            position = int(line[1])
            rs = line[3]
            WT = line[4]
            MUT = line[6]
            position = position - start + 1 # 0-based
            # if the position is recognized as a SNP, ignore
            if position == position_of_interest:
                continue
            else:
                SNP_avoid_range[position] = rs, (WT + "/" + MUT)
    SNPs.close()

####################################################################################################
#################################   Check for secondary structures   ###############################
####################################################################################################
if checks == "yes" or checks =="YES" or checks == "str" or checks == "STR":
    # open the secondary structure file
    sec_str_wt = open(args.u[0])
    sec_str_mut = open(args.U[0])
    # loop all the secondary structures that where found
    for line in sec_str_wt:
        # if the line is not empty
        if line != "":
            line = line.split("\t")
            deltaG_wt = line[2].replace("[", "").replace("]", "")
            sec_str_avoid_range_wt=line[3]
    sec_str_wt.close()
    for line in sec_str_mut:
        # if the line is not empty
        if line != "":
            line = line.split("\t")
            deltaG_mut = line[2].replace("[", "").replace("]", "")
            sec_str_avoid_range_mut=line[3]
    sec_str_mut.close()

####################################################################################################
#################################   Add amplicon and postions found  ###############################
####################################################################################################

# function to find the amplicon and positions
def get_amplicon(template, forward_primer, reverse_primer):
    # Find the start position of the forward primer
    forward_start = template.find(forward_primer)
    if forward_start == -1:
        raise ValueError("Forward primer not found in template sequence")
    
    # Find the start position of the reverse primer (on the reverse complement strand)
    reverse_complement = str(Seq(reverse_primer).reverse_complement())
    reverse_start = template.find(reverse_complement)
    reverse_end = reverse_start + len(reverse_complement)
    if reverse_start == -1:
        raise ValueError("Reverse primer not found in template sequence")
    
    # Extract the amplicon
    amplicon = template[forward_start:reverse_start + len(reverse_complement)]
    return amplicon, forward_start, reverse_end

# function to check the primers
def check_primer_positions(forward_start, reverse_end, SNP_avoid_range, sec_str_avoid_range):
    # check the forward primer
    SNPs_FWD = {}
    sec_FWD = []
    for position,snip in SNP_avoid_range.items():
        if position >= forward_start and position <= (forward_start + len(forward)):
            SNPs_FWD[position-position_of_interest] = snip
    for position in sec_str_avoid_range.replace("[", "").replace("]", "").split(","):
        position = int(position)
        if position >= forward_start and position <= (forward_start + len(forward)):
            sec_FWD.append(position-position_of_interest)
    if not SNPs_FWD:
        SNPs_FWD = "0 found"
    if sec_FWD == []:
        sec_FWD = "0 predicted"
    # check the reverse primer
    SNPs_REV = {}
    sec_REV = []
    for position,snip in SNP_avoid_range.items():
        if position >= reverse_end+1 - len(reverse) and position <= reverse_end:
            SNPs_REV[position-reverse_end] = snip
    for position in sec_str_avoid_range.replace("[", "").replace("]", "").split(","):
        position = int(position)
        if position >= reverse_end+1 - len(reverse) and position <= reverse_end:
            sec_REV.append(position-reverse_end)
    if not SNPs_REV:
        SNPs_REV = "0 found"
    if sec_REV == []:
        sec_REV = "0 predicted"
    return SNPs_FWD, sec_FWD, SNPs_REV, sec_REV

# open the table file
with open(file_output, 'r') as table:
    lines = table.readlines()
with open(file_output, 'w') as output:
    # Header line
# Header line
    output.write("Name\tSpecific_primer\tMatch_Tm\tSingle_MM_Tm\tDouble_MM_Tm\tMM_delta\tGC%\tLenght\tCommon_primer\tMatch_Tm_common\tGC%_common\tLength_common\tFWD_validation\tREV_validation\tSNPs_FWD\tSec_str_FWD\tSNPs_REV\tSec_str_REV\tDeltaG\tAmplicon\tAmp_length\n")    
    for line in lines:
        line = line.rstrip().split("\t")
        # check FWD WT 
        if "_F_WT" in line[0]:
            # get the amplicon and primer positions
            changed_L, changed_R = template_adaptation(left, right ,line[0])
            template = changed_L + wt + right
            forward = line[1]
            reverse = line[8]
            amplicon, forward_start, reverse_end = get_amplicon(template, forward, reverse) #0-based
            # cross-reference the positions
            SNPs_FWD, sec_FWD, SNPs_REV, sec_REV = check_primer_positions(forward_start, reverse_end, SNP_avoid_range, sec_str_avoid_range_wt)
            # write the line to the output file
            output.write("\t".join(line) + "\t" + str(SNPs_FWD) + "\t" + str(sec_FWD) + "\t" + str(SNPs_REV)+ "\t" + str(sec_REV) + "\t" + str(deltaG_wt) + "\t" + str(amplicon) + "\t" + str(len(amplicon)) + "\n")

        # check REV WT
        if "_R_WT" in line[0]:
            # get the amplicon and primer positions
            changed_L, changed_R = template_adaptation(left, right ,line[0])
            template = left + wt + changed_R
            forward = line[8]
            reverse = line[1]
            amplicon, forward_start, reverse_end = get_amplicon(template, forward, reverse) #0-based
            # cross-reference the positions
            SNPs_FWD, sec_FWD, SNPs_REV, sec_REV = check_primer_positions(forward_start, reverse_end, SNP_avoid_range, sec_str_avoid_range_wt)
            # write the line to the output file
            output.write("\t".join(line) + "\t" + str(SNPs_FWD) + "\t" + str(sec_FWD) + "\t" + str(SNPs_REV)+ "\t" + str(sec_REV) + "\t" + str(deltaG_wt) + "\t" + str(amplicon) + "\t" + str(len(amplicon)) + "\n")
        # check FWD MUT
        if "_F_MUT" in line[0]:
            # get the amplicon and primer positions
            changed_L, changed_R = template_adaptation(left, right ,line[0])
            template_MUT = changed_L + m + right
            forward = line[1]
            reverse = line[8]
            amplicon, forward_start, reverse_end = get_amplicon(template_MUT, forward, reverse) #0-based
            # cross-reference the positions
            SNPs_FWD, sec_FWD, SNPs_REV, sec_REV = check_primer_positions(forward_start, reverse_end, SNP_avoid_range, sec_str_avoid_range_mut)
            # write the line to the output file
            output.write("\t".join(line) + "\t" + str(SNPs_FWD) + "\t" + str(sec_FWD) + "\t" + str(SNPs_REV)+ "\t" + str(sec_REV) + "\t" + str(deltaG_mut) + "\t" + str(amplicon) + "\t" + str(len(amplicon)) + "\n")
        # check REV MUT
        if "R_MUT" in line[0]:
            # get the amplicon and primer positions
            changed_L, changed_R = template_adaptation(left, right ,line[0])
            template_MUT = left + m + changed_R
            forward = line[8]
            reverse = line[1]
            amplicon, forward_start, reverse_end = get_amplicon(template_MUT, forward, reverse) #0-based
            # cross-reference the positions
            SNPs_FWD, sec_FWD, SNPs_REV, sec_REV = check_primer_positions(forward_start, reverse_end, SNP_avoid_range, sec_str_avoid_range_mut)
            # write the line to the output file
            output.write("\t".join(line) + "\t" + str(SNPs_FWD) + "\t" + str(sec_FWD) + "\t" + str(SNPs_REV)+ "\t" + str(sec_REV) + "\t" + str(deltaG_mut)  + "\t" + str(amplicon) + "\t" + str(len(amplicon)) + "\n")
table.close()
output.close()

####################################################################################################
#######################################   clean-up the table   #####################################
####################################################################################################

df = pd.read_csv(file_output, sep='\t')
new_order = ["Name","Specific_primer","Match_Tm","Single_MM_Tm","Double_MM_Tm","MM_delta","GC%","Lenght","Common_primer","Match_Tm_common","GC%_common","Length_common","SNPs_FWD","Sec_str_FWD","SNPs_REV","Sec_str_REV","DeltaG","FWD_validation","REV_validation","Amplicon","Amp_length"]
df = df[new_order]
name = ID + "_primer_validated.tsv"
df.to_csv(name, sep='\t', index=False)

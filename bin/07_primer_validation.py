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

parser = argparse.ArgumentParser(description='give arguments to main primer validation script')
parser.add_argument("-F", nargs=1, required=True, help="file with the forward common primers output")
parser.add_argument("-R", nargs=1, required=True, help="file with the reverse common primers output")
parser.add_argument("-T", nargs=1, required=True, help="file with the template sequence; input file")
parser.add_argument("-s", nargs=1, required=True, help="file with the SNP information")
parser.add_argument("-u", nargs=1, required=True, help="file with the secondary structure information")
parser.add_argument("-o", nargs=1, required=True, help="tsv file with the primer information")
parser.add_argument("-p", nargs=1, required=True, help="file with the primer3 settings")
parser.add_argument("-q", nargs=1, required=True, help="upfront checks: yes, snp, str, no")
args = parser.parse_args()

# Get the file names from the arguments
file_forward = args.F[0]
file_reverse = args.R[0]
file_template = args.T[0]
file_snps = args.s[0]
file_structures = args.u[0]
file_output = args.o[0]
primer3_settings= args.p[0]
# get the ID from the file
ID = file_template.split("_")[0].replace(".txt", "")
ID = ID.split("/")[-1]
# upfront check
checks = args.q[0]
SNP_avoid_range = {}
sec_str_avoid_range = []
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
    # The other primers are alternatives and will be written to the alternatives
    elif num != 0:
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
    # The other primers are alternatives and will be written to the alternatives
    elif num != 0:
        alternatives.write(ID + "_common_R_" + str(num) + "\t" + primer + "\t" + tm + "\t" + gc + "\n")
# close the files
Forward.close()
alternatives.close()
Reverse.close()

####################################################################################################
######################################   Validate primers   ########################################
####################################################################################################

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
    for line in lines:
        # Forward primers
        if "_F_WT" in line:
            template = left + wt + right
            # validate the primers
            line_list = line.split("\t")
            name = line_list[0]
            ## needs to be an exact match
            lenght = len(line_list[1])
            temp = left + wt
            forward = temp[-int(lenght):]
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
                line = line.rstrip() + "\t" + common_REV + "\t" + common_TM_REV + "\t" + common_GC_REV + "\t" + str(len(common_REV)) + "\t" + "primer3 failed to validate primers" + "\t" + "primer3 failed to validate primers"+ "\n"
                output.write(line)
        if "_F_MUT" in line:
            template_MUT = left + m + right
            # validate the primers
            line_list = line.split("\t")
            name = line_list[0]
            ## needs to be an exact match
            template_MUT = left + m + right
            lenght = len(line_list[1])
            temp = left + m
            forward = temp[-int(lenght):]
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
                line = line.rstrip() + "\t" + common_REV + "\t" + common_TM_REV + "\t" + common_GC_REV + "\t" + str(len(common_REV)) + "\t" + "primer3 failed to validate primers" + "\t" + "primer3 failed to validate primers"+ "\n"
                output.write(line)
        # Reverse primers
        if "_R_WT" in line:
            template = left + wt + right
            # validate the primers
            line_list = line.split("\t")
            name = line_list[0]
            ## needs to be an exact match
            lenght = len(line_list[1])
            temp = str(Seq(wt + right).reverse_complement())
            reverse = temp[-lenght:]
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
                line = line.rstrip() + "\t" + common_FWD + "\t" + common_TM + "\t" + common_GC + "\t" + str(len(common_FWD)) + "\t" + "primer3 failed to validate primers" + "\t" + "primer3 failed to validate primers"+ "\n"
                output.write(line)
        if "_R_MUT" in line:
            template_MUT = left + m + right
            # validate the primers
            line_list = line.split("\t")
            name = line_list[0]
            ## needs to be an exact match
            lenght = len(line_list[1])
            temp = str(Seq(m + right).reverse_complement())
            reverse = temp[-lenght:]
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
    sec_str = open(args.u[0])
    # loop all the secondary structures that where found
    for line in sec_str:
        # if the line is not empty
        if line != "":
            line = line.split("\t")
            structure = line[1]
        # make a list with the positions to avoid
        i = 0
        for char in list(structure):
            if char == "(":
                sec_str_avoid_range.append(i)
            elif char == ")":
                sec_str_avoid_range.append(i)
            i += 1
    sec_str.close()

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
            SNPs_FWD[position] = snip
    for position in sec_str_avoid_range:
        if position >= forward_start and position <= (forward_start + len(forward)):
            sec_FWD.append(position)
    if not SNPs_FWD:
        SNPs_FWD = "0 found"
    if sec_FWD == []:
        sec_FWD = "0 predicted"
    # check the reverse primer
    positions_REV = range(reverse_end+1 - len(reverse), reverse_end)
    SNPs_REV = {}
    sec_REV = []
    for position,snip in SNP_avoid_range.items():
        if position >= reverse_end+1 - len(reverse) and position <= reverse_end:
            SNPs_REV[position] = snip
    for position in sec_str_avoid_range:
        if position >= reverse_end+1 - len(reverse) and position <= reverse_end:
            sec_REV.append(position)
    if not SNPs_REV:
        SNPs_REV = "0 found"
    if sec_REV == []:
        sec_REV = "0 predicted"
    return SNPs_FWD, sec_FWD, SNPs_REV, sec_REV

# open the table file
with open(file_output, 'r') as table:
    lines = table.readlines()
with open(file_output, 'w') as output:
    for line in lines:
        line = line.rstrip().split("\t")
        # check FWD WT 
        if "F_WT" in line[0]:
            # get the amplicon and primer positions
            template = left + wt + right
            forward = left + wt
            forward = forward[-len(line[1]):]
            reverse = line[8]
            amplicon, forward_start, reverse_end = get_amplicon(template, forward, reverse) #0-based
            # cross-reference the positions
            SNPs_FWD, sec_FWD, SNPs_REV, sec_REV = check_primer_positions(forward_start, reverse_end, SNP_avoid_range, sec_str_avoid_range)
            # write the line to the output file
            output.write("\t".join(line) + "\t" + str(SNPs_FWD) + "\t" + str(sec_FWD) + "\t" + str(SNPs_REV)+ "\t" + str(sec_REV) + "\t" + str(amplicon) + "\n")
        # check REV WT
        if "R_WT" in line[0]:
            # get the amplicon and primer positions
            template = left + wt + right
            forward = line[8]
            reverse = wt + right
            reverse = str(Seq(reverse[0:len(line[1])]).reverse_complement())
            amplicon, forward_start, reverse_end = get_amplicon(template, forward, reverse) #0-based
            # cross-reference the positions
            SNPs_FWD, sec_FWD, SNPs_REV, sec_REV = check_primer_positions(forward_start, reverse_end, SNP_avoid_range, sec_str_avoid_range)
            # write the line to the output file
            output.write("\t".join(line) + "\t" + str(SNPs_FWD) + "\t" + str(sec_FWD) + "\t" + str(SNPs_REV)+ "\t" + str(sec_REV) + "\t" + str(amplicon) + "\n")
        # check FWD MUT
        if "F_MUT" in line[0]:
            # get the amplicon and primer positions
            template_MUT = left + m + right
            forward = left + m
            forward = forward[-len(line[1]):]
            reverse = line[8]
            amplicon, forward_start, reverse_end = get_amplicon(template_MUT, forward, reverse) #0-based
            # cross-reference the positions
            SNPs_FWD, sec_FWD, SNPs_REV, sec_REV = check_primer_positions(forward_start, reverse_end, SNP_avoid_range, sec_str_avoid_range)
            # write the line to the output file
            output.write("\t".join(line) + "\t" + str(SNPs_FWD) + "\t" + str(sec_FWD) + "\t" + str(SNPs_REV)+ "\t" + str(sec_REV) + "\t" + str(amplicon) + "\n")
        # check REV MUT
        if "R_MUT" in line[0]:
            # get the amplicon and primer positions
            template_MUT = left + m + right
            forward = line[8]
            reverse = m + right
            reverse = str(Seq(reverse[0:len(line[1])]).reverse_complement())
            amplicon, forward_start, reverse_end = get_amplicon(template_MUT, forward, reverse) #0-based
            # cross-reference the positions
            SNPs_FWD, sec_FWD, SNPs_REV, sec_REV = check_primer_positions(forward_start, reverse_end, SNP_avoid_range, sec_str_avoid_range)
            # write the line to the output file
            output.write("\t".join(line) + "\t" + str(SNPs_FWD) + "\t" + str(sec_FWD) + "\t" + str(SNPs_REV)+ "\t" + str(sec_REV) + "\t" + str(amplicon) + "\n")
        

table.close()
output.close()
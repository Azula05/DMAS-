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
# get the template sequence
for line in template_file:
    if line != "":
        left = line.split("[")[0]
        wt = line.split("[")[1].split("]")[0].split("/")[0]
        m = line.split("]")[0].split("/")[1]
        right = line.split("]")[1].split("\t")[0]

# append the common priemers to the line and write the new line to the file
with open(file_output,'w') as output:
    for line in lines:
        # Forward primers
        if "_F_WT" in line:
            # validate the primers
            line_list = line.split("\t")
            name = line_list[0]
            ## needs to be an exact match
            lenght = len(line_list[1])
            temp = left + wt + right
            forward = temp[-int(lenght):]
            # validation
            command = "SEQUENCE_ID=" + name + "\nSEQUENCE_TEMPLATE=" + template + "\nSEQUENCE_PRIMER=" + forward + "\nSEQUENCE_PRIMER_REVCOMP=" + common_REV + "\nPRIMER_TASK=check_primers" + "\nPRIMER_EXPLAIN_FLAG=1 \nPRIMER_MIN_TM=45 \nPRIMER_MAX_TM=65 \nPRIMER_MIN_SIZE=16 \nPRIMER_MAX_SIZE=36 \n="
            try:
                process = os.popen("echo \"" + command + "\" | primer3_core -p3_settings_file=" + primer3_settings)            
                validation = process.read()
                process.close()  # Ensure proper resource management
                # Extract stdout from the command
                arguments = validation.split("\n")
                left_validation = arguments[10].split("=")[1]
                right_validation = arguments[11].split("=")[1]
                # change the line
                line = line.rstrip() + "\t" + common_REV + "\t" + common_TM_REV + "\t" + common_GC_REV + "\t" + left_validation + "\t" + right_validation + "\n"
                output.write(line)
            except:
                line = line.rstrip() + "\t" + common_REV + "\t" + common_TM_REV + "\t" + common_GC_REV + "\t" + "primer3 failed to validate primers" + "\t" + "primer3 failed to validate primers"+ "\n"
                output.write(line)
        if "_F_MUT" in line:
            # validate the primers
            line_list = line.split("\t")
            name = line_list[0]
            ## needs to be an exact match
            template_MUT = left + m + right
            lenght = len(line_list[1])
            temp = left + m + right
            forward = temp[-int(lenght):]
            # validation
            command = "SEQUENCE_ID=" + name + "\nSEQUENCE_TEMPLATE=" + template_MUT + "\nSEQUENCE_PRIMER=" + forward + "\nSEQUENCE_PRIMER_REVCOMP=" + common_REV + "\nPRIMER_TASK=check_primers" + "\nPRIMER_EXPLAIN_FLAG=1 \nPRIMER_MIN_TM=45 \nPRIMER_MAX_TM=65 \nPRIMER_MIN_SIZE=16 \nPRIMER_MAX_SIZE=36 \n="
            try:
                process = os.popen("echo \"" + command + "\" | primer3_core -p3_settings_file=" + primer3_settings)            
                validation = process.read()
                process.close()  # Ensure proper resource management
                # Extract stdout from the command
                arguments = validation.split("\n")
                left_validation = arguments[10].split("=")[1]
                right_validation = arguments[11].split("=")[1]
                # change the line
                line = line.rstrip() + "\t" + common_REV + "\t" + common_TM_REV + "\t" + common_GC_REV + "\t" + left_validation + "\t" + right_validation + "\n"
                output.write(line)
            except:
                line = line.rstrip() + "\t" + common_REV + "\t" + common_TM_REV + "\t" + common_GC_REV + "\t" + "primer3 failed to validate primers" + "\t" + "primer3 failed to validate primers"+ "\n"
                output.write(line)
        # Reverse primers
        if "_R_WT" in line:
            # validate the primers
            line_list = line.split("\t")
            name = line_list[0]
            ## needs to be an exact match
            lenght = len(line_list[1])
            temp = str(Seq(left + wt + right).reverse_complement())
            reverse = temp[0:lenght+1]
            # validation
            command = "SEQUENCE_ID=" + name + "\nSEQUENCE_TEMPLATE=" + template + "\nSEQUENCE_PRIMER=" + common_FWD + "\nSEQUENCE_PRIMER_REVCOMP=" + reverse + "\nPRIMER_TASK=check_primers" + "\nPRIMER_EXPLAIN_FLAG=1 \nPRIMER_MIN_TM=45 \nPRIMER_MAX_TM=65 \nPRIMER_MIN_SIZE=16 \nPRIMER_MAX_SIZE=36 \n="
            try:
                process = os.popen("echo \"" + command + "\" | primer3_core -p3_settings_file=" + primer3_settings)            
                validation = process.read()
                process.close()  # Ensure proper resource management
                # Extract stdout from the command
                arguments = validation.split("\n")
                left_validation = arguments[10].split("=")[1]
                right_validation = arguments[11].split("=")[1]
                # change the line
                line = line.rstrip() + "\t" + common_FWD + "\t" + common_TM + "\t" + common_GC + "\t" + left_validation + "\t" + right_validation + "\n"
                output.write(line)
            except:
                line = line.rstrip() + "\t" + common_REV + "\t" + common_TM_REV + "\t" + common_GC_REV + "\t" + "primer3 failed to validate primers" + "\t" + "primer3 failed to validate primers"+ "\n"
                output.write(line)
        if "_R_MUT" in line:
            # validate the primers
            line_list = line.split("\t")
            name = line_list[0]
            ## needs to be an exact match
            lenght = len(line_list[1])
            temp = str(Seq(left + m + right).reverse_complement())
            reverse = temp[0:lenght+1]
            # validation
            command = "SEQUENCE_ID=" + name + "\nSEQUENCE_TEMPLATE=" + template_MUT + "\nSEQUENCE_PRIMER=" + common_FWD + "\nSEQUENCE_PRIMER_REVCOMP=" + reverse + "\nPRIMER_TASK=check_primers" + "\nPRIMER_EXPLAIN_FLAG=1 \nPRIMER_MIN_TM=45 \nPRIMER_MAX_TM=65 \nPRIMER_MIN_SIZE=16 \nPRIMER_MAX_SIZE=36 \n="
            try:
                process = os.popen("echo \"" + command + "\" | primer3_core -p3_settings_file=" + primer3_settings)            
                validation = process.read()
                process.close()  # Ensure proper resource management
                # Extract stdout from the command
                arguments = validation.split("\n")
                left_validation = arguments[10].split("=")[1]
                right_validation = arguments[11].split("=")[1]
                # change the line
                line = line.rstrip() + "\t" + common_FWD + "\t" + common_TM + "\t" + common_GC + "\t" + left_validation + "\t" + right_validation + "\n"
                output.write(line)
            except:
                line = line.rstrip() + "\t" + common_REV + "\t" + common_TM_REV + "\t" + common_GC_REV + "\t" + "primer3 failed to validate primers" + "\t" + "primer3 failed to validate primers"+ "\n"
                output.write(line)
output.close()



####################################################################################################
#######################################   Check for SNPs   #########################################
####################################################################################################

####################################################################################################
#################################   Check for secondary structures   ###############################
####################################################################################################

####################################################################################################
######################################   Add amplicon   ############################################
####################################################################################################
table.close()
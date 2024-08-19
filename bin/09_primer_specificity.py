#!/usr/bin/env python3
"""
This script will check the primers and their specificity to the reference genome.
It will output the results to a new column in the output table.
"""
import argparse
import os
import re

# get the arguments
parser = argparse.ArgumentParser(description="give arguments to filter script")
parser.add_argument('-i', nargs=1, required=True, help="input tsv file with the primers table")
# parser.add_argument('-t', nargs=1, required=True, help="number of threads to use") cannot be used in nextflow with parallel processes (will ask for too many threads)
parser.add_argument('-b', nargs=1, required=True, help="path to the bowtie2 index")
parser.add_argument('-t', nargs=1, required=True, help="Number of threads to use for bowtie2")
parser.add_argument('-s', nargs=1, required=True, help="specificity filter: off, strict, loose")

args = parser.parse_args()
primers_file = args.i[0]
# threads = args.t[0]
threads = args.t[0]
# index
bowtie2_index = args.b[0]
# specificity filter
specificity = args.s[0]


################################################################################################
###################################   Bowtie2   ################################################
################################################################################################
# get the information from the table
# # Open the table and save the contents
with open(primers_file,'r') as table:
	lines = table.readlines()
table.close()
# rewrite the table with the new columns
output_file = open(primers_file, 'w')
# wtite the header
output_file.write("Name\tSpecific_primer\tMatch_Tm\tSingle_MM_Tm\tDouble_MM_Tm\tMM_delta\tGC%\tLenght\tCommon_primer\tMatch_Tm_common\tGC%_common\tLength_common\tSNPs_FWD\tSec_str_FWD\tSNPs_REV\tSec_str_REV\tDeltaG_template\tFWD_validation\tREV_validation\tAmplicon\tAmp_length\tPredicted_structure\tAmplicon_delta_G\tForward_specificity\tReverse_specificity\n")    

if specificity != "off":
	# loop the lines
	line_nr = 0
	for line in lines:
		line = line.strip().split('\t')
		name = line[0]
		## Ignore the header
		if line_nr == 0:
			line_nr += 1
			continue
		## get the primer to check for specificity
		elif "_F_" in line[0]:
			forward = line[1]
			reverse = line[8]
		elif "_R_" in line[0]:
			forward = line[8]
			reverse = line[1]
		## create the input line to check for specificity
		input_forward= ">" + name + "_forward" + "\n" + forward + "\n"
		## run bowtie2
		command = "echo \"" + input_forward +"\" | bowtie2 --no-hd --xeq --no-sq -X 1000 -N 1 --mp 1,1 --quiet -x " + bowtie2_index +" -X1000 --very-sensitive -f - --threads " + threads
		## capture the output
		process = os.popen(command)
		output = process.read()
		## Ensure proper resource management
		process.close()  
		## Extract stdout from the command
		arguments = output.split("\n")
		## Loop all output lines

		# FORWARD
		matches_forward = []
		for argument in arguments:
			## Skip empty lines
			if argument != "":
				## 
				argument = argument.split("\t")
				name = argument[0]
				chromosome = argument[2]
				start = argument[3]
				mismatch = argument[5]
				pattern = re.compile('NM:i:\d+')
				for item in argument:
					if pattern.match(item):
						nm = item.split(":")[2]
				matches_forward.append(chromosome + ":" + start + ":" + mismatch + ":" + nm)
		## create the input line to check for specificity
		input_reverse= ">" + name + "_reverse" + "\n" + reverse
		## run bowtie2
		command = "echo \"" + input_reverse +"\" | bowtie2 --no-hd --xeq --no-sq -X 1000 -N 1 --mp 1,1 --quiet -x " + bowtie2_index +" -X1000 --very-sensitive -f - --threads " + threads
		## capture the output
		process = os.popen(command)
		output = process.read()
		## Ensure proper resource management
		process.close()  
		## Extract stdout from the command
		arguments = output.split("\n")
		## Loop all output lines

		# REVERSE
		matches_reverse = []
		for argument in arguments:
			## Skip empty lines
			if argument != "":
				## 
				argument = argument.split("\t")
				name = argument[0]
				chromosome = argument[2]
				start = argument[3]
				mismatch = argument[5]
				pattern = re.compile('NM:i:\d+')
				for item in argument:
					if pattern.match(item):
						nm = item.split(":")[2]
				matches_reverse.append(chromosome + ":" + start + ":" + mismatch+ ":" + nm)
		## write the output to the table
		output_file.write("\t".join(line) + "\t" + str(matches_forward) + "\t" + str(matches_reverse) + "\n")
line_nr = 0
if specificity == "off":
	for line in lines:
		if line_nr == 0:
			line_nr += 1
			continue
		else:
			output_file.write(line.rstrip() + "\t" + "NA" + "\t" + "NA"+ "\n")
output_file.close()

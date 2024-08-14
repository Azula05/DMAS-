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

args = parser.parse_args()
primers_file = args.i[0]
# threads = args.t[0]
bowtie2_index = args.b[0]
bowtie1_index = args.b[0]
columba_index = args.b[0]

####################################################################################################
#################################    A. specificity filter  ########################################
####################################################################################################

# rules from https://academic.oup.com/clinchem/article/59/10/1470/5622018?login=true fig 6 are used as filter criteria
# The criteria are:

# Strict: 
# - no mismatches in either primer = off-target and is discarted
# - primers with at least 4 mismatches for a single primer = no off-target and is kept
# - primers with a total of at least 5 mismatches between both primers = no off-target and is kept

# Loose:
# - primers with at least 3 mismatches for a single primer = no off-target and is kept
# - primers with a total of at least 4 mismatches between both primers = no off-target and is kept

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

# open the primer file
primers_file = open(primers_file,'r')

# create a dictionary with the primers
primers={}
line_nr = 0
for line in primers_file:
	line = line.strip().split('\t')
	## Ignore the header
	if line_nr == 0:
		line_nr += 1
		continue
	## first one should be the forward primer and the second one the reverse primer
	elif "_F_" in line[0]:
		primers[line[0]] = line[1], line[8]
	elif "_R_" in line[0]:
		primers[line[0]] = line[8], line[1]

primers_file.close()




temp = open("specificity.txt", "w")




################################################################################################
###################################   Bowtie2   ################################################
################################################################################################

# Go through the dictionary and check the specificity of the primers
specificity = {}
i = 0
for name, primer_pair in primers.items():
	input_line = ">" + name + "_forward" + "\n" + primer_pair[0] + "\n" + ">" + name + "_reverse" + "\n" + primer_pair[1]
# echo -e ">primer1_forward\nACTGACTGACTGACTG\n>primer1_reverse\nTGACTGACTGACTGACT" | bowtie2 --threads 2 --no-hd --xeq --no-sq --quiet -x ./Assets/GRCh38/Index_bowtie/GRCh38_noalt_as -f - --very-sensitive -N 1
	# bowtie2 --no-hd --xeq --no-sq -X 1000 -N 1 --mp 1,1 --quiet -x " + bowtie2_index +" -X1000 --very-sensitive -f - --threads 3
	# bowtie --tryhard -X1000 -v3 --quiet -x " + bowtie1_index + " --quiet --threads 3 -f -
	command = "echo \"" + input_line +"\" | bowtie2 --no-hd --xeq --no-sq -X 1000 -N 1 --mp 1,1 --quiet -x " + bowtie2_index +" -X1000 --very-sensitive -f - --threads 3"
	print(i)
	temp.write(str(i)+"\n")
	process = os.popen(command)
	output = process.read()
	process.close()  # Ensure proper resource management
	# Extract stdout from the command
	arguments = output.split("\n")
	temp.write(str(arguments)+"\n")
	# USE RE WILL NOT WORK LIKE THIS
	#specificity[name] = arguments[0].split("\t")[2], arguments[0].split("\t")[3], arguments[0].split("\t")[5], arguments[0].split("\t")[17].split(":")[2] ,arguments[1].split("\t")[2], arguments[1].split("\t")[3], arguments[1].split("\t")[5], arguments[1].split("\t")[17].split(":")[2]
	#print(specificity)
	#if i == 6:
		#exit("break")
	i += 1

## ADD CPUS AAND MAKE  SURE IIIIIT IS NOT PARAALLEL AAND USES AALL


"""
################################################################################################
###################################   Bowtie1   ################################################
################################################################################################
specificity = {}
i = 0
for name, primer_pair in primers.items():
	input_line = ">" + name + "_forward" + "\n" + primer_pair[0] + "\n" + ">" + name + "_reverse" + "\n" + primer_pair[1]
# echo -e ">primer1_forward\nACTGACTGACTGACTG\n>primer1_reverse\nTGACTGACTGACTGACT" | bowtie2 --threads 2 --no-hd --xeq --no-sq --quiet -x ./Assets/GRCh38/Index_bowtie/GRCh38_noalt_as -f - --very-sensitive -N 1
	# "echo \"" + input_line +"\" | bowtie2 --no-hd --xeq --no-sq --quiet -x " + bowtie2_index +" -f - --very-sensitive -N 1"
	command = "echo \"" + input_line +"\" | bowtie --tryhard -X1000 -v3 --quiet -x " + bowtie1_index + " --quiet --threads 3 -f -"
	print(i)
	temp.write(str(i))
	process = os.popen(command)
	output = process.read()
	process.close()  # Ensure proper resource management
	# Extract stdout from the command
	arguments = output.split("\n")
	temp.write(str(arguments)+"\n")
	# USE RE WILL NOT WORK LIKE THIS
	#specificity[name] = arguments[0].split("\t")[2], arguments[0].split("\t")[3], arguments[0].split("\t")[5], arguments[0].split("\t")[17].split(":")[2] ,arguments[1].split("\t")[2], arguments[1].split("\t")[3], arguments[1].split("\t")[5], arguments[1].split("\t")[17].split(":")[2]
	#print(specificity)
	i += 1
temp.close()
"""
"""
################################################################################################
###################################   columba   ################################################
################################################################################################

specificity = {}
i = 0
for name, primer_pair in primers.items():
	input_line = ">" + name + "_forward" + "\n" + primer_pair[0] + "\n" + ">" + name + "_reverse" + "\n" + primer_pair[1]
# echo -e ">primer1_forward\nACTGACTGACTGACTG\n>primer1_reverse\nTGACTGACTGACTGACT" | bowtie2 --threads 2 --no-hd --xeq --no-sq --quiet -x ./Assets/GRCh38/Index_bowtie/GRCh38_noalt_as -f - --very-sensitive -N 1
	# "echo \"" + input_line +"\" | bowtie2 --no-hd --xeq --no-sq --quiet -x " + bowtie2_index +" -f - --very-sensitive -N 1"
	command = "echo \"" + input_line +"\" | columba primer_search --primers - --targets genome.fasta"
	print(i)
	process = os.popen(command)
	output = process.read()
	process.close()  # Ensure proper resource management
	# Extract stdout from the command
	arguments = output.split("\n")
	print(arguments)
	# USE RE WILL NOT WORK LIKE THIS
	#specificity[name] = arguments[0].split("\t")[2], arguments[0].split("\t")[3], arguments[0].split("\t")[5], arguments[0].split("\t")[17].split(":")[2] ,arguments[1].split("\t")[2], arguments[1].split("\t")[3], arguments[1].split("\t")[5], arguments[1].split("\t")[17].split(":")[2]
	#print(specificity)
	if i == 6:
		exit("break")
	i += 1
"""
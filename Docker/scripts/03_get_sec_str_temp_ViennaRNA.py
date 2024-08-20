#!/usr/bin/env python3
"""
This script is part of the main DMAS pipeline and is used to get the secondary structure of the template sequences using ViennaRNA.
It requires the user to provide the template sequence -t
"""
# ARGUMENT HANDLING
import argparse
import os

parser = argparse.ArgumentParser(description="give arguments to main DMAS script")
parser.add_argument("-t", nargs=1, required=True, help="sequence")
args = parser.parse_args()

#######################################################################################################################
################   GET SECONDARY STRUCTURE OF THE TEMPLATE SEQUENCES USING VIENNARNA   ################################
#######################################################################################################################

# SEQUENCE ID
seq_ID = os.path.splitext(args.t[0])[0]
seq_ID = seq_ID.split("/")[-1]

# SEQUENCE
seq = open(args.t[0]).readline().rstrip().split('\t')[0]

# retrieve the wild type of the seq (whole sequence with [N/])
wild = seq.split("[")[0] + seq.split("[")[1].split("/")[0] + seq.split("]")[1]
# retrieve the mutant type of the seq (whole sequence with [/N])
mutant = seq.split("[")[0] + seq.split("/")[1].split("]")[0] + seq.split("]")[1]

############################################################################################
##################################   Wild type   ###########################################
############################################################################################
command = "echo " + wild +"| RNAfold -p --noconv -T 60.0 -P DNA --salt 0.05"
process = os.popen(command)
output = process.read()
process.close()  # Ensure proper resource management
arguments = output.split("\n")

## structure
WT_structure = arguments[3].split(" ")[0]
## delta_g
WT_delta_g = arguments[2].split(" ",1)[1]

############################################################################################
##################################    Mutation   ###########################################
############################################################################################
command = "echo " + mutant +"| RNAfold -p --noconv -T 60.0 -P DNA --salt 0.05"
process = os.popen(command)
output = process.read()
process.close()  # Ensure proper resource management
arguments = output.split("\n")

## structure
MUT_structure = arguments[3].split(" ")[0]
## delta_g
MUT_delta_g = arguments[2].split(" ",1)[1]

############################################################################################
#################################    Output files   ########################################
############################################################################################

# Wild type
output_wt = open("sec-str_" + seq_ID + '_wt' + ".txt", "a")
output_wt.write(seq_ID + '_wt' + '\t' + WT_structure + "\t" + str(WT_delta_g) + "\t")
# Which positions have secondary structures?
str_not_ok = []
# if seq has secondary structures ("(" or ")"), retrieve the positions and add it to a list
for index, j in enumerate(WT_structure):
    if j == "(" or j == ")":
        str_not_ok.append(index) # 0-based
output_wt.write(str(str_not_ok) + "\n")
output_wt.close()

## Mutant
output_mut = open("sec-str_" + seq_ID + '_mut'+ ".txt", "a")
output_mut.write( seq_ID + '_var-mut'+ '\t' + MUT_structure + "\t" + str(MUT_delta_g) + "\t")
# Which positions have secondary structures?
str_not_ok = []
# if seq has secondary structures ("(" or ")"), retrieve the positions and add it to a list
for index, j in enumerate(MUT_structure):
    if j == "(" or j == ")":
        str_not_ok.append(index) # 0-based
output_mut.write(str(str_not_ok) + "\n")
output_mut.close()


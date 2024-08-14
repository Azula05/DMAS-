#!/usr/bin/env python3
"""
This script will predict the secondary structure of the amplicons.
"""
import argparse
import os

# get the arguments
parser = argparse.ArgumentParser(description="give arguments to main DMAS script")
parser.add_argument('-i', nargs=1, required=True, help="input tsv file with the primers table")
args = parser.parse_args()

input_file = args.i[0]

# Get the ID
id = input_file.split('_')[0]
id = id.split("/")[-1]
id = id .replace('.tsv', '')

# Get the amplicons
## Get the information that is in the file (back-up)
with open(input_file, 'r') as input:
    lines = input.readlines()
input.close()
## Rewrite the file
with open(input_file, 'w') as output:
    line_nr = 0
    for line in lines:
        # Write the header
        if line_nr == 0:
            output.write("Name\tSpecific_primer\tMatch_Tm\tSingle_MM_Tm\tDouble_MM_Tm\tMM_delta\tGC%\tLenght\tCommon_primer\tMatch_Tm_common\tGC%_common\tLength_common\tFWD_validation\tREV_validation\tSNPs_FWD\tSec_str_FWD\tSNPs_REV\tSec_str_REV\tDeltaG\tAmplicon\tAmp_length\tPredicted_structure\tAmplicon_delta_G\n")    
        # Write the rest of the lines
        else:
            ### Get the amplicon
            line = line.rstrip().split('\t')
            try:
                amplicon = line[19]
                ### Run viennaRNA
                command = "echo " + amplicon +"| RNAfold -p --noconv -T 60.0 -P DNA --salt 0.05"
                process = os.popen(command)
                output2 = process.read()
                process.close()
                arguments = output2.split("\n")
                ### Get the MFE proxy structures
                structure = str(arguments[3].split(" ")[0])
                ### Get the delta G
                delta_g = arguments[2].split(" ",1)[1]
                delta_g = delta_g.replace("[","").replace("]","").strip()
                output.write("\t".join(line) + "\t" + structure + "\t" + delta_g + "\n")
            except:
                structure = "NA"
                delta_g = "NA"
                output.write("\t".join(line) + "\t" + structure + "\t" + delta_g + "\n")
        line_nr += 1


output.close()

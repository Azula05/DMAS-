#!/usr/bin/env python3

from nupack import *
import argparse
import os

parser = argparse.ArgumentParser(description="give arguments to main DMAS script")
parser.add_argument("-t", nargs=1, required=True, help="sequence")
args = parser.parse_args()


# get seq ID from file name
seq_ID = os.path.splitext(args.t[0])[0]
# get sequence from input file
seq = open(args.t[0]).readline().rstrip().split('\t')[0]


# retrieve the wild type of the seq (whole sequence with [N/])
wild = seq.split("[")[0] + seq.split("[")[1].split("/")[0] + seq.split("]")[1]
# retrieve the mutant type of the seq (whole sequence with [/N])
mutant = seq.split("[")[0] + seq.split("/")[1].split("]")[0] + seq.split("]")[1]

# NUPACK model initiation
mdl = Model(material = "DNA", celsius = 60, sodium = 0.05, magnesium = 0.003)

# compute MFE proxy structures
# wt
my_mfe = mfe(strands = wild, model = mdl)
delta_g = float(my_mfe[0].energy)
structure = str(my_mfe[0].structure)
# mut
my_mfe_mut = mfe(strands = mutant, model = mdl)
delta_g_mut = float(my_mfe[0].energy)
structure_mut = str(my_mfe[0].structure)

# which sequence has the worst secondary structures? wt or mut?
# write structure and deltag of worst sequence to output file
if structure == structure_mut:
    output = open("sec-str_" + seq_ID + '_var-wtmt' + ".txt", "a")
    output.write(seq_ID + '_var-wtmt' + '\t' + structure + "\t" + str(delta_g) + "\t")
elif (structure.count("(") + structure.count(")")) > (structure_mut.count("(") + structure_mut.count(")")):
    output = open("sec-str_" + seq_ID + '_var-wt'+ ".txt", "a")
    output.write( seq_ID + '_var-wt'+ '\t' + structure + "\t" + str(delta_g) + "\t")
elif (structure.count("(") + structure.count(")")) < (structure_mut.count("(") + structure_mut.count(")")):
    output = open("sec-str_" + seq_ID + '_var-mt'+ ".txt", "a")
    structure = structure_mut
    delta_g = delta_g_mut
    output.write(seq_ID + '_var-mt'+ '\t' + structure + "\t" + str(delta_g) + "\t")

str_not_ok = []
# if seq has secondary structures ("(" or ")"), retrieve the positions and add it to a list
for index, j in enumerate(structure):
    if j == "(" or j == ")":
        str_not_ok.append(index)

# write list with positions of secondary structures to output file
# generate the output file

output.write(str(str_not_ok) + "\n")

output.close()

#!/usr/bin/env python3

from nupack import *
import argparse
import os

parser = argparse.ArgumentParser(description="give arguments to main DMAS script")
args = parser.parse_args()

file_list = os.listdir()
file_list_in = []
print(file_list)

for file_name in file_list:
    if file_name[0:22] == "input-sec-str-amplicon":
        file_list_in.append(file_name)

print(file_list_in)

for file_name in file_list_in:
    seq_ID = file_name.split('_')[1]
    templ_ID = file_name.replace('.txt', '').split('_')[2]

    # check if there where any primers
    if len((file_name).split('_')) > 3:
        primer_ID = file_name.split('_')[3]

        output = open("out-sec-str-amp_" + seq_ID + '_' + templ_ID + '_' + primer_ID + ".txt", "a")

        amp_file = open(file_name).readline()
        tmp, seq_ID, templ_ID, primer_ID, amplicon = amp_file.rstrip().split("_")
            
        # NUPACK model initiation
        mdl = Model(material = "DNA", celsius = 60, sodium = 0.05, magnesium = 0.003)

        # compute MFE proxy structures
        my_mfe = mfe(strands = amplicon, model = mdl)
        delta_g = float(my_mfe[0].energy)
        structure = str(my_mfe[0].structure)

        # make string for regions to avoid
        str_not_ok = []
        for index, i in enumerate(structure):
            if i == "(" or i == ")":
                str_not_ok.append(index)

        output.write(seq_ID + "_" + templ_ID + "_" + primer_ID + "\t" + str(delta_g) + "\t" + str(str_not_ok) + "\t" + structure + "\n")

    else:
        output = open("out-sec-str-amp_" + seq_ID + '_' + templ_ID  + ".txt", "a")
        output.write(seq_ID + '_' + templ_ID + '\t' + "no_primers_found\n")



    output.close()

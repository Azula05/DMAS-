#!/usr/bin/env python3

import argparse
import os

parser = argparse.ArgumentParser(description="give arguments to main DMAS script")
args = parser.parse_args()

# check files in directory
file_list = os.listdir()
file_list_in = []

# get all files with sec-str-amplicon in the name
for file_name in file_list:
    if file_name[0:22] == "input-sec-str-amplicon":
        file_list_in.append(file_name)

# get sec-str-amplicon for each file
for file_name in file_list_in:
    seq_ID = file_name.split('_')[1]
    templ_ID = file_name.replace('.txt', '').split('_')[2]
    
    # check if there where any primers
    if len((file_name).split('_')) > 3:
        primer_ID = file_name.split('_')[3]
        primer_ID = primer_ID.replace('.txt', '')
        output = open("out-sec-str-amp_" + seq_ID + '_' + templ_ID + '_' + primer_ID + ".txt", "a")

        amp_file = open(file_name).readline()
        tmp, seq_ID, templ_ID, primer_ID, amplicon = amp_file.rstrip().split("_")
        
        # Do prediction of secondary structure
        command = "echo " + amplicon +"| RNAfold -p --noconv -T 60.0 -P DNA --salt 0.05"
        process = os.popen(command)
        output2 = process.read()
        process.close()  # Ensure proper resource management
        arguments = output2.split("\n")

        # compute MFE proxy structures
        delta_g = arguments[2].split(" ",1)[1]
        structure = str(arguments[3].split(" ")[0])

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

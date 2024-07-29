#!/usr/bin/env python3

import argparse
import re
import os


parser = argparse.ArgumentParser(description="give arguments to filter script")
parser.add_argument("-i", nargs = 1, required=True, help="file with sequence info")
parser.add_argument("-p", nargs = 1, required=True, help="file with all designed primers")
parser.add_argument("-s", nargs = 1, required=True, help="file with primer specificity")
parser.add_argument("-S", nargs=1, required=True, help="bed file with SNPs") # nothing happens with this info
parser.add_argument("-t", nargs=1, required=True, help="file with secondary structures") # nothing happes with this info

args = parser.parse_args()

# this script runs on a sequence level and gathers and adds filter info for all primers for that specific sequence


#### Get all possible template IDs form sec_str.txt file and make a list for which no primers where found

sec_str_amp = open(args.t[0])

no_design_list = set()
design_list = set()
amp_fold_dict = {}

for line in sec_str_amp:
    if len(line.split('\t')) == 2:
        no_design_list.add(line.split('\t')[0])
    else:
        seq_ID, templ_ID, primer_ID = line.split('\t')[0].split('_')
        design_list.add(seq_ID + "_" + templ_ID)
        amp_fold_dict[line.split('\t')[0]] = line.rstrip().split('\t')[1:]


# file with all designed primers for this sequence
all_primers = open(args.p[0])           
seq_ID = os.path.splitext(args.i[0])[0]

# make a new file for the filtered primers
filtered_primers = open("filtered-primers_" + seq_ID +  ".txt", "w")

#### Get specificity info
spec_file = open(args.s[0]).readlines()      # file with primer specificity info
spec_list = []
for line in spec_file:
    spec_list.append(line.rstrip())

# make counters for log file
passed = 0
failed_spec = 0
failed_str_amp = 0
design = 1
total_primers = 0
# warning_snp = 0
# warning_sec_str = 0

# get mutant variation

seq = open(args.i[0]).readline().split('\t')[0]

c = 0
for nt in seq:
    if nt == "/":
        mutant = seq[c+1]
        wt_mt = seq[c-1:c+2]
    c += 1

# if no primers could be designed
if os.path.getsize(args.p[0]) == 0:
    design = 0

# rev comp dict
comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}


# go through all the primers and filter them

primer_found = 0

for primer_pair in all_primers:
    total_primers += 1
    # get all the info
    seq_ID, templ_ID, primer_ID, f_seq, r_seq, f_pos, f_len, r_pos, r_len, f_tm, r_tm, f_GC, r_GC, amplicon  = primer_pair.split()

    f_pos = int(f_pos)
    f_len = int(f_len)
    r_pos = int(r_pos)
    r_len = int(r_len)

    complete_ID = seq_ID + '_' + templ_ID + '_' + primer_ID
    
    filter_str = ""
    warning_str = ""

    # off targets
    if complete_ID in spec_list:
        filter_str = filter_str + "FAIL_specificity_"

    # get amplicon folding info into dict
    amp_fold_avoid = {}
    amp_fold_avoid_pos = {}
    deltag, str_pos, str_pos_nr = amp_fold_dict[complete_ID]
    amp_fold_avoid[primer_ID] = float(deltag)
    amp_fold_avoid_pos[primer_ID] = str_pos


    fold_amp_avoid = amp_fold_avoid_pos[primer_ID]

    if fold_amp_avoid != '[]':
        fold_amp_avoid = fold_amp_avoid.replace('[', "").replace(']', '').split(', ')
        fold_amp_avoid = [ int(x) for x in fold_amp_avoid ]

    if amp_fold_avoid[primer_ID] < -15:
        filter_str = filter_str + 'FAIL_fold_amplicon_'
    elif (amp_fold_avoid[primer_ID] < -5) & any(x in range(f_len + 1) for x in fold_amp_avoid):
        filter_str = filter_str + 'FAIL_fold_amplicon_'
    elif (amp_fold_avoid[primer_ID] < -5) & any(x in range(r_pos - r_len - f_pos, r_pos - f_pos + 1) for x in fold_amp_avoid):
        filter_str = filter_str + 'FAIL_fold_amplicon_'
    


    # if all tests succeded => primer pair passed filters
    if filter_str == "":
        filter_str = "PASS_"
        primer_found = 1

    # gather info to add to log file
    if filter_str == "PASS_":
        passed += 1
    if "FAIL_specificity" in filter_str:
        failed_spec += 1
    if "FAIL_fold_amplicon" in filter_str:
        failed_str_amp += 1

    # generate DMAS primer

    if templ_ID[6:9] == 'fwd':
        mutant_seq = f_seq[:-1] + mutant
    elif templ_ID[6:9] == 'rev':
        mutant_seq = r_seq[:-1] + comp[mutant]

    filtered_primers.write(seq_ID + "\t" + templ_ID + "\t" + primer_ID + "\t" + seq + '\t' + wt_mt + '\t' + f_seq + "\t" + r_seq + "\t" + mutant_seq + '\t' + str(f_pos) + "\t" + str(f_len) + "\t" + str(r_pos) + "\t" + str(r_len) + "\t" + f_tm + "\t" + r_tm + "\t" + f_GC + "\t" + r_GC + "\t" + amplicon + "\t" + filter_str[0:len(filter_str)-1] + "\n")

all_primers.close()

# what if no primers for none of the options could be designed
if design == 0:
    filtered_primers.write(seq_ID + "\tno primers could be designed for this sequence\n")

elif primer_found == 0:
    filtered_primers.write(seq_ID + "\tno primer pair passed all filters for this sequence\n")


filtered_primers.close()


# print log file
log_file = open("log-file_" + seq_ID + ".txt", "a")

log_file.write(
    seq_ID + '\t' + 
    str(design) + "\t" + str(primer_found) + "\t" + str(total_primers) + "\t" + 
    str(passed) + '\t' + str(failed_spec) + "\t" + str(failed_str_amp) + '\n' )
log_file.close()

#!/usr/bin/env python3

import os
import argparse

parser = argparse.ArgumentParser(description="give arguments to filter script")
parser.add_argument("-f", nargs = 1, required=True, help="file with filtered primers")

#file_list = os.listdir()
#file_list_primers = [f for f in file_list if f[:16] == "filtered-primers"] # only keep the files we need

primer_found_list = set()

output_file = open("dmas_primers.txt", "w")
output_file.write("seq_ID\ttempl_ID\tprimer_ID\tinput_seq\twt_mt\tf_seq\tr_seq\tdmas_seq\tf_pos\tf_length\tr_pos\tr_length\tf_tm\tr_tm\tf_GC\tr_GC\tamplicon\n")

args = parser.parse_args()

primer_file = open(args.f[0])

for primer_line in primer_file:
    if len(primer_line.split('\t')) == 2:
        output_file.write(primer_line)
        
    else:
        seq_ID, templ_ID, primer_ID = primer_line.split("\t")[:3]

        primer_filter = primer_line.rstrip().split('\t')[-1]

        if seq_ID not in primer_found_list: # if no primer has been selected for this input sequence
            if primer_filter == "PASS":

                output_file.write("\t".join(primer_line.split('\t')[:-1]) + '\n')
                primer_found_list.add(seq_ID)

primer_file.close()


output_file.close()


#summary file

#file_list_logs = [f for f in file_list if f[:16] == "filtered-primers"] # only keep the files we need







# from itertools import groupby



# if primers[0] != "": # not sure why this if statement is here => for when there is no selected primer pair?!
#     dict = {}
#     c = 0
#     sequence = []
#     for i in primers:
#         # the first line contains the header, so it will be skipped
#         if c == 0:
#             c += 1
#         else:
#             if i[0] not in sequence:
#                 sequence.append(i[0])
            
#             temp_ID, ID, chrom, start, end, primer_ID, f_seq, r_seq, f_pos, f_len, r_pos, r_len, f_tm, r_tm, f_GC, r_GC, a, warning_str, fail_str, sec_str_amp, f_MM, r_MM, amp_snp, sec_str_temp  = i.rstrip().split("\t")
#             # initialize the score (number of mismatches/snps)
#             score = 0
#             # add amount of forward & reverse primer mismatches to score
#             score += int(f_MM)
#             score += int(r_MM)
#             # when amplicon contains no snps where primer anneals, the list is empty
#             if (amp_snp.replace("[", "").replace("]", "").split(",")) == [""]:
#                 # don't do anything
#                 score += 0
#             # when the list is not empty, add number of snps to score
#             else:
#                 score += len(amp_snp.replace("[", "").replace("]", "").split(","))
#             # when template contains no secondary structures where primer anneals, the list is empty
#             if (sec_str_temp.replace("[", "").replace("]", "").split(",")) == [""]:
#                 # don't do anything
#                 score += 0
#             # when the list is not empty, add number of secondary structure positions to score
#             else:
#                 score += len(sec_str_temp.replace("[", "").replace("]", "").split(","))
#             # for all primers that passed filtering steps:
#             if "PASS" in fail_str:
#                 # write a dictionary with all information as key and the score as the value
#                 dict[temp_ID + "\t" + ID + "\t" + chrom + "\t" + start + "\t" + end + 
#                      "\t" + primer_ID + "\t" + f_seq + "\t" + r_seq + "\t" + f_pos + 
#                      "\t" + f_len + "\t" + r_pos + "\t" + r_len + "\t" + f_tm + 
#                      "\t" + r_tm + "\t" + f_GC + "\t" + r_GC + "\t" + a + "\t" + warning_str + 
#                      "\t" + sec_str_amp + "\t" + f_MM + "\t" + r_MM + "\t" + amp_snp + 
#                      "\t" + sec_str_temp + "\n"] = score
#         c += 1
#     sequence.sort()
#     if c == 2:
#         for i in sequence:
#             temp_ID, ID, chrom, start, end, primer_ID, f_seq, r_seq, f_pos, f_len, r_pos, r_len, f_tm, r_tm, f_GC, r_GC, a, warning_str, sec_str_amp, f_MM, r_MM, amp_snp, sec_str_temp  = list(dict.keys())[0].rstrip().split("\t")
#             if i == temp_ID:
#                 # for forward allele-specific primers
#                 if int(ID) < 3:
#                     nature = "FWD"
#                     MM_pos = int(ID) + 2    # mismatch position
#                 # for reverse allele-specific primers
#                 else:
#                     nature = "REV"
#                     MM_pos = int(ID) - 1        # mismatch position
#                 warning = ""
#                 if "WARNING_snp" in warning_str:
#                     warning = warning + "Common snp at AS primer position"
#                 if "WARNING_sec_str" in warning_str:
#                     # if there is already a warning, provide a comma to separate them
#                     if warning != "":
#                         warning = warning + ", Secondary structure at AS primer position"
#                     else:
#                         warning = warning + "Secondary structure at AS primer position"
#                 output.write(temp_ID + "\t" + f_seq + "\t" + r_seq + "\t" + nature + "\t" + str(MM_pos) + "\t" + f_MM + "\t" + r_MM + "\t" + warning + "\n")
        
#     else:
#         # sort dictionary on score (ascending)
#         score_sorted = {k: v for k, v in sorted(dict.items(), key=lambda item: item[1])}

#         # get the first value of the dictionary
#         dict_values = score_sorted.values()
#         iterator = iter(dict_values)
#         value = next(iterator)

#         # create a subset from the dictionary for all items that have the same value as the first one
#         subset = {k: v for k, v in dict.items() if v == value}
#         # if there is only one item with that score
#         for i in sequence:
#             if len(subset) == 1:
#                 temp_ID, ID, chrom, start, end, primer_ID, f_seq, r_seq, f_pos, f_len, r_pos, r_len, f_tm, r_tm, f_GC, r_GC, a, warning_str, sec_str_amp, f_MM, r_MM, amp_snp, sec_str_temp  = list(subset.keys())[0].rstrip().split("\t")
#                 if i == temp_ID:
#                     # for forward allele-specific primers
#                     if int(ID) < 3:
#                             nature = "FWD"
#                             MM_pos = int(ID) + 2    # mismatch position
#                     # for reverse allele-specific primers
#                     else:
#                         nature = "REV"
#                         MM_pos = int(ID) - 1        # mismatch position
#                     warning = ""
#                     if "WARNING_snp" in warning_str:
#                         warning = warning + "Common snp at AS primer position"
#                     if "WARNING_sec_str" in warning_str:
#                         # if there is already a warning, provide a comma to separate them
#                         if warning != "":
#                             warning = warning + ", Secondary structure at AS primer position"
#                         else:
#                             warning = warning + "Secondary structure at AS primer position"
#                     output.write(temp_ID + "\t" + f_seq + "\t" + r_seq + "\t" + nature + "\t" + str(MM_pos) + "\t" + f_MM + "\t" + r_MM + "\t" + warning + "\n")
#             # if there is more then 1 item with the same score as the first one
#             else:
#                 dict = {}
#                 for j in subset:
#                     temp_ID, ID, chrom, start, end, primer_ID, f_seq, r_seq, f_pos, f_len, r_pos, r_len, f_tm, r_tm, f_GC, r_GC, a, warning_str, sec_str_amp, f_MM, r_MM, amp_snp, sec_str_temp  = j.rstrip().split("\t")
#                     # write a dictionary with all information as key and the amplicon length as the value
#                     dict[temp_ID + "\t" + ID + "\t" + chrom + "\t" + start + "\t" + end + 
#                         "\t" + primer_ID + "\t" + f_seq + "\t" + r_seq + "\t" + f_pos + 
#                         "\t" + f_len + "\t" + r_pos + "\t" + r_len + "\t" + f_tm + 
#                         "\t" + r_tm + "\t" + f_GC + "\t" + r_GC + "\t" + a + "\t" + warning_str + 
#                         "\t" + sec_str_amp + "\t" + f_MM + "\t" + r_MM + "\t" + amp_snp + 
#                         "\t" + sec_str_temp + "\n"] = len(a)
                
#                 # sort dictionary on amplicon length (ascending)
#                 amp_sorted = {k: v for k, v in sorted(dict.items(), key=lambda item: item[1])}
#                 MM_pos_list = []
#                 for j in list(amp_sorted.keys()):
#                     temp_ID, ID, chrom, start, end, primer_ID, f_seq, r_seq, f_pos, f_len, r_pos, r_len, f_tm, r_tm, f_GC, r_GC, a, warning_str, sec_str_amp, f_MM, r_MM, amp_snp, sec_str_temp  = j.rstrip().split("\t")
#                     if i == temp_ID:
#                         # if a primer with that certain artificial mismatch position is not already written to the output file
#                         if ID not in MM_pos_list:
#                             MM_pos_list.append(ID)
#                             # for forward allele-specific primers
#                             if int(ID) < 3:
#                                     nature = "FWD"
#                                     MM_pos = int(ID) + 2    # mismatch position
#                             # for reverse allele-specific primers
#                             else:
#                                 nature = "REV"
#                                 MM_pos = int(ID) - 1        # mismatch position
#                             warning = ""
#                             if "WARNING_snp" in warning_str:
#                                 warning = warning + "Common snp at AS primer position"
#                             if "WARNING_sec_str" in warning_str:
#                                 # if there is already a warning, provide a comma to separate them
#                                 if warning != "":
#                                     warning = warning + ", Secondary structure at AS primer position"
#                                 else:
#                                     warning = warning + "Secondary structure at AS primer position"
#                             output.write(temp_ID + "\t" + f_seq + "\t" + r_seq + "\t" + nature + "\t" + str(MM_pos) + "\t" + f_MM + "\t" + r_MM + "\t" + warning + "\n")

# output.close()
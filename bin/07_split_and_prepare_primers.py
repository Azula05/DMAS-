#!/usr/bin/env python3

import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to main snp script')
args = parser.parse_args()


file_list = os.listdir()
file_list_in = []

for file_name in file_list:
	if file_name[0:13] == 'input-primer3':
		file_list_in.append(file_name)


for file_name in file_list_in:
	os.system("/bin/primer3-2.6.1/src/primer3_core --output=out-" + file_name + " --p3_settings_file=${settings_ch} " + file_name)



for file_name in file_list_in:
	primers_in = open('out-' + file_name)

	# convert primer3 output content to dictionary
	primers = {}

	for line in primers_in:
		key, value = line.rstrip().split("=")
		primers[key] = value

	primers_in.close()

	# get relevant output from created dictionary
	temp_ID, seq_ID = (primers["SEQUENCE_ID"]).split("_")
	template = primers["SEQUENCE_TEMPLATE"]
	nr_p_out = int(primers["PRIMER_LEFT_NUM_RETURNED"])

	# read general info into file
	info_keys = ("SEQUENCE_ID", "SEQUENCE_TEMPLATE", "SEQUENCE_TARGET")
	general_info = open("primer-design-info_" + temp_ID  + '_' + seq_ID + ".txt", "a")

	for info in primers:
		if "_NUM_" in info or "_EXPLAIN" in info or any(x in info for x in info_keys):
			general_info.write(info + '=' + str(primers[info]) +'\n')

	general_info.close()

	# make file for bowtie
	bowtie_file = open("input-primer-spec_" + temp_ID  + '_' + seq_ID  + ".txt", "a")

	# # make general file with list primers
	all_primers_file = open("primer-list_" + temp_ID + '_' + seq_ID + ".txt", 'w')
	all_primers_dict = {}

	if nr_p_out > 0:
		for primer_index in range(nr_p_out):


			FWD = primers[("PRIMER_LEFT_" + str(primer_index) + "_SEQUENCE")]
			FWD_qual = len(FWD) * "I"
			REV = primers[("PRIMER_RIGHT_" + str(primer_index) + "_SEQUENCE")]
			REV_qual = len(REV) * "I"

			PRIMER_LEFT_TM = primers[("PRIMER_LEFT_" + str(primer_index) + "_TM")]
			PRIMER_RIGHT_TM = primers[("PRIMER_RIGHT_" + str(primer_index) + "_TM")]
			PRIMER_LEFT_GC_PERCENT = primers[("PRIMER_LEFT_" + str(primer_index) + "_GC_PERCENT")]
			PRIMER_RIGHT_GC_PERCENT = primers[("PRIMER_RIGHT_" + str(primer_index) + "_GC_PERCENT")]

			# bowtie input file
			# write FWD + REV
			bowtie_file.write(temp_ID  + '_' + seq_ID + "_primer-" + str(primer_index) + "_FWD_REV" + "\t")
			bowtie_file.write(FWD + "\t" + FWD_qual + "\t" + REV + "\t" + REV_qual + "\n")

			# write  REV + FWD
			bowtie_file.write(temp_ID  + '_' + seq_ID + "_primer-" + str(primer_index) + "_REV_FWD" + "\t")
			bowtie_file.write(REV + "\t" + REV_qual + "\t" + FWD + "\t" + FWD_qual + "\n")


			# write FWD + FWD
			bowtie_file.write(temp_ID  + '_' + seq_ID + "_primer-" + str(primer_index) + "_FWD_FWD" + "\t")
			bowtie_file.write(FWD + "\t" + FWD_qual + "\t" + FWD + "\t" + FWD_qual + "\n")

			# write REV + REV
			bowtie_file.write(temp_ID  + '_' + seq_ID + "_primer-" + str(primer_index) + "_REV_REV" + "\t")
			bowtie_file.write(REV + "\t" + REV_qual + "\t" + REV + "\t" + REV_qual + "\n")


			# NUPACK files
			amplicon_file = open("input-sec-str-amplicon_" + temp_ID + "_" + seq_ID + '_primer-' + str(primer_index) + ".txt", "w")

			# write lines to NUPACK input file
			FWD_pos, FWD_len = primers['PRIMER_LEFT_'+ str(primer_index)].split(",")
			REV_pos, REV_len = primers['PRIMER_RIGHT_'+ str(primer_index)].split(",")
			amplicon = template[int(FWD_pos):int(REV_pos) + 1]									 # get the amplicon using info from primer3 output file
			amplicon_file.write("amplicon_" + temp_ID + "_" + seq_ID + "_primer-" + str(primer_index) + "_" + amplicon + "\n")
			
			amplicon_file.close()


			# general primer file (for filtering), first put in dict, will be sorted (see below)
			all_primers_dict[temp_ID + "\t" + seq_ID + '\tprimer-' + str(primer_index) + '\t' + FWD + '\t' + REV + '\t' + 
			FWD_pos + '\t' + FWD_len + '\t' + REV_pos +'\t' + REV_len + '\t' + PRIMER_LEFT_TM + '\t' + PRIMER_RIGHT_TM + '\t' + 
			PRIMER_LEFT_GC_PERCENT + '\t' + PRIMER_RIGHT_GC_PERCENT + '\t' + amplicon + '\n'] = len(amplicon)

	else:
		amplicon_file = open("input-sec-str-amplicon_" + temp_ID + "_" + seq_ID + ".txt", "w")
		amplicon_file.write("no primers found")
		amplicon_file.close()


	# sort primers according to amp size (smallest is best) and then print to all_primers file
	all_primers_sorted = {k: v for k, v in sorted(all_primers_dict.items(), key=lambda item: item[1])}

	for primer in all_primers_sorted:
		all_primers_file.write(primer)

	bowtie_file.close()
	all_primers_file.close()
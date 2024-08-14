#!/usr/bin/env python3

import argparse
import re
import os


parser = argparse.ArgumentParser(description="give arguments to filter script")
parser.add_argument('-f', nargs=1, required=True, help="stringency of filtering primer specificity")

args = parser.parse_args()
spec_filter = args.f[0]

# this script summarises the SAM output file to give a PASS or FAIL for each primer pair

spec_file = open('out-primer-spec_all.sam')

all_lines = []

for line in spec_file:
	all_lines.append(line)

spec_file.close()


fail_spec_file = open('fail-spec_all.txt', 'w')


# rules from https://academic.oup.com/clinchem/article/59/10/1470/5622018?login=true fig 6
for line_nr in range(0, len(all_lines) - 1, 2):
	seq_ID = all_lines[line_nr].split()[0].split('_')[0]
	templ_ID = all_lines[line_nr].split()[0].split('_')[1]
	primer_ID = all_lines[line_nr].split()[0].split('_')[2]

	fwd_spec = all_lines[line_nr].split()
	rev_spec = all_lines[line_nr+1].split()

	# get info about that seq ID
	seq_file = open(seq_ID + ".txt").readline().strip()
	seq, seq_pos, SNP_pos = seq_file.split('\t')
	chrom, start, end = seq_pos.replace(":", "-").split("-")

	# check if matches are region of interest
	if fwd_spec[2] == chrom and int(fwd_spec[3]) - 1 <= int(end) and int(fwd_spec[3]) - 1 >= int(start):
		continue										# when the chromosome in the specificity file is the same one as in the primer file and its forward position falls between sequence positions: continue
	if rev_spec[2] == chrom and int(rev_spec[3]) - 1 <= int(end) and int(rev_spec[3]) - 1 >= int(start):
		continue										# when the chromosome in the specificity file is the same one as in the primer file and its reverse position falls between sequence positions: continue

	else:
		fwd_MM, rev_MM = 0, 0

		if len(fwd_spec) > 7:
			fwd_MM = fwd_spec[7].count('>')

		if len(rev_spec) > 7:
			rev_MM = rev_spec[7].count('>')

		if fwd_MM > 0 or rev_MM > 0:
			if spec_filter == 'strict':
				if fwd_MM + rev_MM < 5:
					fail_spec_file.write(seq_ID + '_' + templ_ID + '_' + primer_ID + '\n')
			if spec_filter == 'loose':
				if (fwd_MM + rev_MM < 3) or (fwd_MM == 1 & rev_MM == 2) or (fwd_MM == 2 & rev_MM == 1):
					fail_spec_file.write(seq_ID + '_' + templ_ID + '_' + primer_ID + '\n')
		else:
			fail_spec_file.write(seq_ID + '_' + templ_ID + '_' + primer_ID + '\n')

fail_spec_file.close()
#!/usr/bin/env python3

# Arguments
import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to main snp script')
parser.add_argument("-a", nargs=1, required=True, help="primer min left 3 prime distance")
parser.add_argument("-b", nargs=1, required=True, help="number of primers to return")
parser.add_argument("-c", nargs=1, required=True, help="min melting temperature")
parser.add_argument("-d", nargs=1, required=True, help="max melting temperature")
parser.add_argument("-e", nargs=1, required=True, help="optimal melting temperature")
parser.add_argument("-f", nargs=1, required=True, help="melting temperature difference between primers")
parser.add_argument("-g", nargs=1, required=True, help="min GC content")
parser.add_argument("-i", nargs=1, required=True, help="max GC content")
parser.add_argument("-j", nargs=1, required=True, help="optimal GC content")
parser.add_argument("-k", nargs=1, required=True, help="min amplicon length")
parser.add_argument("-l", nargs=1, required=True, help="max amplicon length")
parser.add_argument("-m", nargs=1, required=True, help="mispriming library")
parser.add_argument("-M", nargs=1, required=True, help="max mispriming library")
parser.add_argument("-q", nargs=1, required=True, help="variable for filtering options")
parser.add_argument("-s", nargs=1, required=True, help="bed file with SNPs")
parser.add_argument("-t", nargs=1, required=True, help="sequence")
parser.add_argument("-u", nargs=1, required=True, help="file with secondary structures")
args = parser.parse_args()

# Determine SNP type
def split_up(sequence):
    """
    Pre-processes input into different sequences and determines SNP type
    Argument: DNA template to design primers for
    Returns: splitted up template, SNP and SNP type
    """
    before = sequence.split("/")[0]     # Normal nucleotide
    after = sequence.split("/")[1].rstrip() # Mutant nucleotide
    before, wt = before.split("[")      # sequence before SNP, wild type nucleotide
    length = len(before)                # = position where allele-specific primer anneals
    m, after = after.split("]")         # mutant nucleotide, sequence after SNP
    if len(wt) == len(m):               # determine the SNP type
        snp_type = "exchange"
    elif len(wt) < len(m):
        snp_type = "insertion"
    elif len(wt) > len(m):
        snp_type = "deletion"
    return before, after, length, wt, m, snp_type

# Generate templates for exchange SNP (WT)
def templ_generation_exchange(before, after, wt):
    """
    Creates 6 templates from pre-processed sequence as input for primer3:
        - 3 possible artificial mismatch positions (2, 3, 4)
        - 2 possible AS primer locations (forward, reverse)
    Args: pre-processed sequences and SNP
    Returns: list with adjusted templates
    """    
    # forward allele-specific primer (mismatch created in sequence before SNP)
    seqs = {}
    # nucleotides 2-4 are changed to create mismatch
    for mm_pos in range(2,5):
        seq = list(before)
        # nucleotide is replaced with C, except when it already contains a C,
        # in that case it is replace with A
        if seq[-mm_pos] == "C" or seq[-mm_pos] == "c":
            seq[-mm_pos] = "A"
        else:
            seq[-mm_pos] = "C"
        seqs["templ-fwd-"+str(mm_pos)] = "".join(seq) + wt + after
    
    # reverse allele-specific primer (mismatch created in sequence after SNP)
    # nucleotides 2-4 are changed to create mismatch
    for mm_pos in range(1,4):
        seq = list(after)
        # nucleotide is replaced with G, except when it already contains a G,
        # in that case it is replace with T
        if seq[mm_pos] == "G" or seq[mm_pos] == "g":
            seq[mm_pos] = "T"
        else:
            seq[mm_pos] = "G"
        seqs["templ-rev-"+str(mm_pos+1)] = before + wt + "".join(seq) # +1 as the REV primer is counted in the other direction (does not start at 0)

    return(seqs)


# get all the info from the sequence input file
seq_info = open(args.t[0]).readline().rstrip()
seq, pos, targeted_snp_pos = seq_info.split('\t')
chrom, start, end = pos.replace(':', '-').split('-')
start = int(start)
end = int(end)
targeted_snp_pos = int(targeted_snp_pos)
seq_length = len(seq)
seq_ID = os.path.splitext(args.t[0])[0]
seq_ID = seq_ID.split("/")[-1]


# get separate parts of input sequence
before, after, length, wt, m, snp_type = split_up(seq)

# check the snp_type of the SNP
if snp_type == "exchange":
    # generate the 6 possible template sequences
    seqs = templ_generation_exchange(before, after, wt)
# if snp_type == "insertion":
#   seqs = templ_generation_insertion()
# if snp_type == "deletion":
#   seqs = templ_generation_deletion()


# Make list from start and end positions, including all positions in between
# This is needed to compare with positions of SNPs and return the index of this list
all_pos = list(range(start, end))


# is filter on?
upfront_filter = args.q[0]
# Adding SNP positions to SEQUENCE_EXCLUDED_REGION

if upfront_filter == "yes" or upfront_filter == "snp":
    SNP_all = []
    SNP_file = open(args.s[0])
    for i in SNP_file:                                      # Add positions of SNPs to list
        if chrom == i.split()[0]:
            SNP_pos = int(i.split()[1]) + 1                     # Determine SNP start position
            SNP_length = int(i.split()[2]) - int(i.split()[1])  # Determine SNP length

            for i in range(SNP_length):                         # Depending on the SNP length, add all SNP positions that cover this length (often only length of 1)
                SNP_all.append(SNP_pos)
                SNP_pos += 1
    SNP_file.close()

# Adding secondary structure positions to SEQUENCE_EXCLUDED_REGION
if upfront_filter == "yes" or upfront_filter == "str":                      # If secondary structure filter is on
    sec_str_file = open(args.u[0]).readline().strip()
    fold_temp_avoid = sec_str_file.split('\t')[3]      

    if fold_temp_avoid != "[]":
        fold_temp_avoid = fold_temp_avoid.replace("[", "").replace("]", "").split(", ")
        fold_temp_avoid = [int(x) for x in fold_temp_avoid] # Get list of all secondary structure positions

### Create primer3 input files
for templ_ID, templ_seq in seqs.items():
    output = open("input-primer3_" + seq_ID + '_' + templ_ID + ".txt", "w")
    output.write("SEQUENCE_ID=" + seq_ID + '_' + templ_ID + "\n")
    output.write("SEQUENCE_TEMPLATE=" + templ_seq + "\n")
    output.write("PRIMER_NUM_RETURN=" + args.b[0] + "\n")
    output.write("PRIMER_MIN_THREE_PRIME_DISTANCE=" + args.a[0] + "\n")
    output.write("PRIMER_PRODUCT_SIZE_RANGE=" + args.k[0] + "-" + args.l[0] + "\n")
    output.write("PRIMER_MIN_TM=" + args.c[0] + "\n")
    output.write("PRIMER_MAX_TM=" + args.d[0] + "\n")
    output.write("PRIMER_OPT_TM=" + args.e[0] + "\n")
    output.write("PRIMER_PAIR_MAX_DIFF_TM=" + args.f[0] + "\n")
    output.write("PRIMER_MIN_GC=" + args.g[0] + "\n")
    output.write("PRIMER_MAX_GC=" + args.i[0] + "\n")
    output.write("PRIMER_OPT_GC_PERCENT=" + args.j[0] + "\n")

    if templ_ID[6:9] == 'fwd':  # fwd templates
        output.write("SEQUENCE_FORCE_LEFT_END=" + str(length) + "\n")
    elif templ_ID[6:9] == 'rev':   # last three templates have a MM in REV primer
        output.write("SEQUENCE_FORCE_RIGHT_END=" + str(length) + "\n")



    # make SEQUENCE_EXCLUDED_REGION, if upfront filter is on
    if upfront_filter != "no":

        # when creating forward AS primers + common reverse primers
        if templ_ID[6:9] == 'fwd':  # fwd templates

            avoid_range = []

            for snp_pos in SNP_all:
                if int(all_pos.index(snp_pos)) <= targeted_snp_pos: # If the SNP position comes before the targeted SNP, don't include it in the primer3 setting
                    continue
                else:                                           # If the SNP position comes after the targeted SNP, include it in the primer3 setting
                    avoid_range.append([all_pos.index(snp_pos), 1])      # Append SNPs in format for primer3 setting

            for str_pos in fold_temp_avoid:
                if str_pos <= targeted_snp_pos:              # If the structure position comes before the targeted SNP, don't include it in the primer3 setting
                    continue
                else:                                       # If the structure position comes after the targeted SNP, include it in the primer3 setting
                    avoid_range.append([str_pos, 1])              # Append structure position in format for primer3 setting


            if avoid_range != []:
                avoid_range_str = "SEQUENCE_EXCLUDED_REGION="
                for i in avoid_range:                                   # Write new settings string
                    avoid_range_str = avoid_range_str + str(i[0]) + "," + str(i[1]) + " "
                output.write(avoid_range_str + '\n')


        # when creating reverse AS primers + common forward primers
        elif templ_ID[6:9] == 'rev':   # last three templates have a MM in REV primer
            avoid_range = []
            
            for i in SNP_all:
                if int(all_pos.index(i)) >= targeted_snp_pos: # If the SNP position comes after the targeted SNP, don't include it in the primer3 setting
                    continue
                else:                                           # If the SNP position comes before the targeted SNP, include it in the primer3 setting
                    avoid_range.append([all_pos.index(i), 1])

            for i in fold_temp_avoid:
                if i >= targeted_snp_pos:              # If the structure position comes after the targeted SNP, don't include it in the primer3 setting
                    continue
                else:                                       # If the structure position comes before the targeted SNP, include it in the primer3 setting
                    avoid_range.append([i, 1])

            if avoid_range != []:
                avoid_range_str = "SEQUENCE_EXCLUDED_REGION="
                for i in avoid_range:                                   # Write new settings string
                    avoid_range_str = avoid_range_str + str(i[0]) + "," + str(i[1]) + " "
            
            output.write(avoid_range_str + '\n')


    # screen primers against library on which they cannot anneal (only if 
    # user specified to do so)
    if args.m[0] != "no":
        print(args.m[0])
        output.write("PRIMER_MISPRIMING_LIBRARY=" + args.m[0] + "\n")
        output.write("PRIMER_MAX_LIBRARY_MISPRIMING=" + args.M[0] + "\n")
    output.write("=")
    output.close()

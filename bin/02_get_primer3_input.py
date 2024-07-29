#!/usr/bin/env python3

import argparse
import os

parser = argparse.ArgumentParser(description="give arguments to main primer DMAS script")
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
parser.add_argument("-t", nargs=1, required=True, help="template")
args = parser.parse_args()


def split_up(template):
    """
    Pre-processes input into different sequences and determines SNP type
    Argument: DNA template to design primers for
    Returns: splitted up template, SNP and SNP type
    """
    before = template.split("/")[0]
    after = template.split("/")[1].rstrip()
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


def templates_exchange(before, after, wt):
    """
    Creates 6 sequences from pre-processed template as input for primer3:
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


# open input sequence
templ = open(args.t[0]).readline()
seq_ID = os.path.splitext(args.t[0])[0]

# get separate parts of input sequence
before, after, length, wt, m, snp_type = split_up(templ)

# check the snp_type of the SNP
if snp_type == "exchange":
    # generate the 6 possible template sequences
    seqs = templates_exchange(before, after, wt)
# if snp_type == "insertion":
#   seqs = templates_insertion()
# if snp_type == "deletion":
#   seqs = templates_deletion()


### Create primer3 input file
for templ_ID, templ_seq in seqs.items():
    output = open("input_primer3_" + templ_ID + ".txt", "w")
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
    # if creating forward AS primers + common reverse primers
    if templ_ID[6:9] == 'fwd':  # first three templates have a MM in FWD primer
        output.write("SEQUENCE_FORCE_LEFT_END=" + str(length) + "\n")
    # if creating reverse AS primers + common forward primers
    elif templ_ID[6:9] == 'rev':   # last three templates have a MM in REV primer
        output.write("SEQUENCE_FORCE_RIGHT_END=" + str(length) + "\n")
    # screen primers against library on which they cannot anneal (only if 
    # user specified to do so)
    if args.m[0] != "no":
        print(args.m[0])
        output.write("PRIMER_MISPRIMING_LIBRARY=" + args.m[0] + "\n")
        output.write("PRIMER_MAX_LIBRARY_MISPRIMING=" + args.M[0] + "\n")
    output.write("=")
    output.close()


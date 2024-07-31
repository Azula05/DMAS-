#!/usr/bin/env python3
"""
This script creates input files for primer3 for all possible templates.
It considers SNPs, secondary structures and mispriming libraries if these options are turned on.
All positions are considered for mismatch except when one is specified.
"""
# ARGUMENT HANDLING
import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to main snp script')
parser.add_argument("-a", nargs=1, required=True, help="primer min left 3 prime distance") # Is dit nodig? wat doet dit?
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
parser.add_argument("-q", nargs=1, required=True, help="variable for filtering options: yes, snp, str, no")
parser.add_argument("-s", nargs=1, required=True, help="bed file with SNPs")
parser.add_argument("-t", nargs=1, required=True, help="seq-file with template")
parser.add_argument("-u", nargs=1, required=True, help="file with secondary structures")
parser.add_argument("-p", nargs=1, required=True, help="Position of the mismatch eighter all, 2, 3 or 4 postions to 5' end")
args = parser.parse_args()

MIN_LEFT_3_PRIME_DISTANCE = args.a[0]
NUMBER_OF_PRIMERS = args.b[0]
MIN_MELTING_TEMPERATURE = args.c[0]
MAX_MELTING_TEMPERATURE = args.d[0]
OPTIMAL_MELTING_TEMPERATURE = args.e[0]
MELTING_TEMPERATURE_DIFFERENCE = args.f[0]
MIN_GC_CONTENT = args.g[0]
MAX_GC_CONTENT = args.i[0]
OPTIMAL_GC_CONTENT = args.j[0]
MIN_AMPLICON_LENGTH = args.k[0]
MAX_AMPLICON_LENGTH = args.l[0]
MISPRIMING_LIBRARY = args.m[0]
MAX_MISPRIMING_LIBRARY = args.M[0]
FILTERING_OPTIONS = args.q[0]
SNP_file = args.s[0]
seq_file = args.t[0]
sec_str_file = args.u[0]
position_mismatch = args.p[0]

############################################################################################################
#################################    FUNCTION: SPLIT SEQUENCE    ###########################################
############################################################################################################
def split_up(sequence):
    """
    Pre-processes input into different sequences parts [WT/MUT] and determines SNP type.
    Argument: DNA template to design primers for
    Returns: splitted up template, SNP and SNP type and position in 1-based format
    """
    sequence = sequence.upper()             # make all uppercase
    before = sequence.split("/")[0]         # 5' of the template
    after = sequence.split("/")[1].rstrip() # 3' of the template
    before, wt = before.split("[")          # sequence before SNP, wild type nucleotide
    length = len(before)                    # = position where allele-specific primer anneals (1-based!)
    m, after = after.split("]")             # mutant nucleotide, sequence after SNP
    # determine the SNP type
    if len(wt) == len(m):               
        snp_type = "exchange"
    elif len(wt) < len(m):
        snp_type = "insertion"
    elif len(wt) > len(m):
        snp_type = "deletion"
    return before, after, length, wt, m, snp_type

############################################################################################################
#################################    FUNCTION: Exchange WT    ##############################################
############################################################################################################
# Generate templates for exchange SNP (WT)
def templ_generation_exchange(before, after, wt, position_mismatch):
    """
    Creates 6 templates from pre-processed sequence as input for primer3:
        - 3 possible artificial mismatch positions (2, 3, 4)
        - 2 AS primer (forward, reverse)
    Args: pre-processed sequences and the SNP
    Returns: list with adjusted templates
    """    
    # A) FORWARD PRIMERS
    # forward allele-specific primer (mismatch created in sequence before SNP to the 5' end)
    seqs = {}
    # all 3 positions are considered (2, 3, 4)
    if position_mismatch == "all":
        #A) FORWARD PRIMERS
        for mm_pos in range(2,5):
            seq = list(before) # Everything that comes before the SNP
            # Replace C (C,A,T)
            # Replace G (G,A,T)
            # Replace T (T,G,C)
            # Replace A (A,G,C)

        # B) REVERSE PRIMERS
        for mm_pos in range(2,5): # wrong in the old script (1,2,3 instead of 2,3,4)
            seq = list(after)
            # nucleotide is replaced with G, except when it already contains a G,
            # in that case it is replace with T
            if seq[mm_pos] == "G" or seq[mm_pos] == "g":
                seq[mm_pos] = "T"
            else:
                seq[mm_pos] = "G"
            seqs["templ-rev-"+str(mm_pos+1)] = before + wt + "".join(seq) # +1 as the REV primer is counted in the other direction (does not start at 0)

    return(seqs)

############################################################################################################
#################################    Test   ############################################
############################################################################################################

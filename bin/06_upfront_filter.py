#!/usr/bin/env python3
"""
This script creates input files for primer3 for all possible templates.
It considers SNPs, secondary structures and mispriming libraries if these options are turned on.
All positions are considered for mismatch except when one is specified.
"""
# ARGUMENT HANDLING
import argparse
from Bio.Seq import Seq

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
###################################    FUNCTION: AS 1MM TEMPLATES    #######################################
############################################################################################################

def AS_1MM_templates(before, after, wt, m, position_mismatch, seq_ID):

    # A) Allele specific primers (1MM)
    # A good length for PCR primers is generally around 18-30 bases

    # Function to get the reverse complement of a sequence
    def reversecomplement(seq):
        my_seq = Seq(seq)
        complement = my_seq.reverse_complement()
        return str(complement)
    
    # A) FORWARD AS PRIMERS (1MM)
    AS_FWD_WT = before + wt
    AS_FWD_MUT = before + m

    # B) REVERSE AS PRIMERS (1MM)
    REV_WT = wt + after
    AS_REV_WT = reversecomplement(REV_WT)
    REV_MUT = m + after
    AS_REV_MUT = reversecomplement(REV_MUT)
    AS_primers = {seq_ID + "_AS_FWD_WT" : AS_FWD_WT, seq_ID + "_AS_FWD_MUT" : AS_FWD_MUT, seq_ID + "_AS_REV_WT" : AS_REV_WT, seq_ID + "_AS_REV_MUT" : AS_REV_MUT}
    return(AS_primers)

############################################################################################################
####################################    FUNCTION: FWD templates    #########################################
############################################################################################################

def FORWARD_TEMPLATES(before, wt, m, position_mismatch, seq_ID):
    # FORWARD PRIMERS
    ART_primers = {}
    if position_mismatch == "all":
        pos = [2,3,4]
    else:
        pos = [position_mismatch]
    for mm_pos in pos:
        seq = list(before)          # Everything that comes before the SNP
        Original = seq[-mm_pos]     # Original base
        # Replace C (C,A,T)
        if Original == "C":
            temp = list(before)
            # Replacement
            temp[-mm_pos] = "A"
            # WT FWD
            WT = "".join(temp) + wt
            ART_primers[seq_ID + "_" + str(mm_pos) + "_A_FWD_WT"] = WT
            # MUT FWD
            MUT = "".join(temp) + m
            ART_primers[seq_ID + "_" + str(mm_pos) + "_A_FWD_MUT"] = MUT
        if Original == "C":
            temp = list(before)
            temp[-mm_pos] = "C"
            WT = "".join(temp) + wt
            ART_primers[seq_ID + "_" + str(mm_pos) + "_C_FWD_WT"] = WT
            # MUT FWD
            MUT = "".join(temp) + m
            ART_primers[seq_ID + "_" + str(mm_pos) + "_C_FWD_MUT"] = MUT
        if Original == "C":
            temp = list(before)
            temp[-mm_pos] = "T"
            WT = "".join(temp) + wt
            ART_primers[seq_ID + "_" + str(mm_pos) + "_T_FWD_WT"] = WT
            # MUT FWD
            MUT = "".join(temp) + m
            ART_primers[seq_ID + "_" + str(mm_pos) + "_T_FWD_MUT"] = MUT         
        # Replace G (G,A,T)
        if Original == "G":
            temp = list(before)
            temp[-mm_pos] = "G"
            WT = "".join(temp) + wt
            ART_primers[seq_ID + "_" + str(mm_pos) + "_G_FWD_WT"] = WT
            # MUT FWD
            MUT = "".join(temp) + m
            ART_primers[seq_ID + "_" + str(mm_pos) + "_G_FWD_MUT"] = MUT
        if Original == "G":
            temp = list(before)
            temp[-mm_pos] = "A"
            WT = "".join(temp) + wt
            ART_primers[seq_ID + "_" + str(mm_pos) + "_A_FWD_WT"] = WT
            # MUT FWD
            MUT = "".join(temp) + m
            ART_primers[seq_ID + "_" + str(mm_pos) + "_A_FWD_MUT"] = MUT
        if Original == "G":
            temp = list(before)
            temp[-mm_pos] = "T"
            WT = "".join(temp) + wt
            ART_primers[seq_ID + "_" + str(mm_pos) + "_T_FWD_WT"] = WT
            # MUT FWD
            MUT = "".join(temp) + m
            ART_primers[seq_ID + "_" + str(mm_pos) + "_T_FWD_MUT"] = MUT   
        # Replace T (T,G,C)
        if Original == "T":
            temp = list(before)
            temp[-mm_pos] = "T"
            WT = "".join(temp) + wt
            ART_primers[seq_ID + "_" + str(mm_pos) + "_T_FWD_WT"] = WT
            # MUT FWD
            MUT = "".join(temp) + m
            ART_primers[seq_ID + "_" + str(mm_pos) + "_T_FWD_MUT"] = MUT
        if Original == "T":
            temp = list(before)
            temp[-mm_pos] = "G"
            WT = "".join(temp) + wt
            ART_primers[seq_ID + "_" + str(mm_pos) + "_G_FWD_WT"] = WT
            # MUT FWD
            MUT = "".join(temp) + m
            ART_primers[seq_ID + "_" + str(mm_pos) + "_G_FWD_MUT"] = MUT
        if Original == "T":
            temp = list(before)
            temp[-mm_pos] = "C"
            WT = "".join(temp) + wt
            ART_primers[seq_ID + "_" + str(mm_pos) + "_C_FWD_WT"] = WT
            # MUT FWD
            MUT = "".join(temp) + m
            ART_primers[seq_ID + "_" + str(mm_pos) + "_C_FWD_MUT"] = MUT      
        # Replace A (A,G,C)
        if Original == "A":
            temp = list(before)
            temp[-mm_pos] = "A"
            WT = "".join(temp) + wt
            ART_primers[seq_ID + "_" + str(mm_pos) + "_A_FWD_WT"] = WT
            # MUT FWD
            MUT = "".join(temp) + m
            ART_primers[seq_ID + "_" + str(mm_pos) + "_A_FWD_MUT"] = MUT
        if Original == "A":
            temp = list(before)
            temp[-mm_pos] = "G"
            WT = "".join(temp) + wt
            ART_primers[seq_ID + "_" + str(mm_pos) + "_G_FWD_WT"] = WT
            # MUT FWD
            MUT = "".join(temp) + m
            ART_primers[seq_ID + "_" + str(mm_pos) + "_G_FWD_MUT"] = MUT
        if Original == "A":
            temp = list(before)
            temp[-mm_pos] = "C"
            WT = "".join(temp) + wt
            ART_primers[seq_ID + "_" + str(mm_pos) + "_C_FWD_WT"] = WT
            # MUT FWD
            MUT = "".join(temp) + m
            ART_primers[seq_ID + "_" + str(mm_pos) + "_C_FWD_MUT"] = MUT
    return ART_primers

############################################################################################################
####################################    FUNCTION: Exchange    ##############################################
############################################################################################################

# Generate templates for exchange SNP
#def templ_generation_exchange(before, after, wt, m, mm_pos, seq_ID):
    """
    Creates templates from pre-processed sequence as input for primer3:
        - 2 AS primer (forward, reverse) => 1MM
        - 3 possible artificial mismatch positions (2, 3, 4) => create long primers for temperature and primer validation
    Args: pre-processed sequences and the SNP
    Returns: list with adjusted templates
    """


############################################################################################################
#################################    Test   ############################################
############################################################################################################

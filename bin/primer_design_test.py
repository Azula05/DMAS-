#!/usr/bin/env python3
from Bio.Seq import Seq
sequence = "TGTTATGAACAGCAAATAAAAGAAACTAAAACGATCCTGAGACTTCCACACTGATGCAATCATTCGTCTGTTTCCCATTCTAAACTGTACCCTGTTACTT[A/C]TCCCCTTCCTATGACATGAACTTAACCATAGAAAAGAAGGGGAAAGAAAACATCAAGCGTCCCATAGACTCACCCTGAAGTTCTCAGGATCCACGTGCAG"

before = "TGTTATGAACAGCAAATAAAAGAAACTAAAACGATCCTGAGACTTCCACACTGATGCAATCATTCGTCTGTTTCCCATTCTAAACTGTACCCTGTTACTT"
after = "TCCCCTTCCTATGACATGAACTTAACCATAGAAAAGAAGGGGAAAGAAAACATCAAGCGTCCCATAGACTCACCCTGAAGTTCTCAGGATCCACGTGCAG"
wt = "A"
m = "C"
position_mismatch = "all"
seq_ID = "seq-0"
"""
Creates templates from pre-processed sequence as input for primer3:
    - 2 AS primer (forward, reverse) => 1MM
    - 3 possible artificial mismatch positions (2, 3, 4) => create long primers for temperature and primer validation
Args: pre-processed sequences and the SNP
Returns: list with adjusted templates
"""

def templ_generation_exchange(before, after, wt, m, position_mismatch, seq_ID):
    """
    Creates templates from pre-processed sequence as input for primer3:
        - 2 AS primer (forward, reverse) => 1MM
        - 3 possible artificial mismatch positions (2, 3, 4) => create long primers for temperature and primer validation
    Args: pre-processed sequences and the SNP
    Returns: list with adjusted templates
    """
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



AS_primers = templ_generation_exchange(before, after, wt, m, position_mismatch, seq_ID)
print(AS_primers)


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

ART_primers_FWD = FORWARD_TEMPLATES(before, wt, m, position_mismatch, seq_ID)

print(ART_primers_FWD)

def REVERSE_TEMPLATES(after, wt, m, position_mismatch, seq_ID):
    ART_primers = {}
    if position_mismatch == "all":
        pos = [2,3,4]
    else:
        pos = [position_mismatch]
    for mm_pos in pos:
        seq = list(before)          # Everything that comes before the SNP
        Original = seq[-mm_pos]     # Original base
        # Replace C (C,A,T)
    return ART_primers
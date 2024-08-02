#!/usr/bin/env python3
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

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
print("AS_primers:")
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

print("ART_primers_FWD:")
print(ART_primers_FWD)

def REVERSE_TEMPLATES(after, wt, m, position_mismatch, seq_ID):
    ART_primers = {}
    def reversecomplement(seq):
        my_seq = Seq(seq)
        complement = my_seq.reverse_complement()
        return str(complement)
    if position_mismatch == "all":
        pos = [2,3,4]
    else:
        pos = [position_mismatch]
    for mm_pos in pos:
            seq = list(after)          # Everything that comes before the SNP
            Original = seq[mm_pos]     # Original base
            # Replace C (C,A,T)
            if Original == "C":
                temp = list(after)
                # WT REV
                WT = wt + "".join(temp)
                WT = reversecomplement(WT)
                WT = list(WT)
                WT[-(mm_pos+1)] = "C"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_C_REV_WT"] = "".join(WT)
                # MUT REV
                MUT = m + "".join(temp)
                MUT = reversecomplement(MUT)
                MUT = list(MUT)
                MUT[-(mm_pos+1)] = "C"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_C_REV_MUT"] = "".join(MUT)
            if Original == "C":
                temp = list(after)
                # WT REV
                WT = wt + "".join(temp)
                WT = reversecomplement(WT)
                WT = list(WT)
                WT[-(mm_pos+1)] = "A"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_A_REV_WT"] = "".join(WT)
                # MUT REV
                MUT = m + "".join(temp)
                MUT = reversecomplement(MUT)
                MUT = list(MUT)
                MUT[-(mm_pos+1)] = "A"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_A_REV_MUT"] = "".join(MUT)
            if Original == "C":
                temp = list(after)
                # WT REV
                WT = wt + "".join(temp)
                WT = reversecomplement(WT)
                WT = list(WT)
                WT[-(mm_pos+1)] = "T"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_T_REV_WT"] = "".join(WT)
                # MUT REV
                MUT = m + "".join(temp)
                MUT = reversecomplement(MUT)
                MUT = list(MUT)
                MUT[-(mm_pos+1)] = "T"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_T_REV_MUT"] = "".join(MUT)
            # Replace G (G,A,T)
            if Original == "G":
                temp = list(after)
                # WT REV
                WT = wt + "".join(temp)
                WT = reversecomplement(WT)
                WT = list(WT)
                WT[-(mm_pos+1)] = "G"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_G_REV_WT"] = "".join(WT)
                # MUT REV
                MUT = m + "".join(temp)
                MUT = reversecomplement(MUT)
                MUT = list(MUT)
                MUT[-(mm_pos+1)] = "G"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_G_REV_MUT"] = "".join(MUT)
            if Original == "G":
                temp = list(after)
                # WT REV
                WT = wt + "".join(temp)
                WT = reversecomplement(WT)
                WT = list(WT)
                WT[-(mm_pos+1)] = "A"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_A_REV_WT"] = "".join(WT)
                # MUT REV
                MUT = m + "".join(temp)
                MUT = reversecomplement(MUT)
                MUT = list(MUT)
                MUT[-(mm_pos+1)] = "A"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_A_REV_MUT"] = "".join(MUT)
            if Original == "G":
                temp = list(after)
                # WT REV
                WT = wt + "".join(temp)
                WT = reversecomplement(WT)
                WT = list(WT)
                WT[-(mm_pos+1)] = "T"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_T_REV_WT"] = "".join(WT)
                # MUT REV
                MUT = m + "".join(temp)
                MUT = reversecomplement(MUT)
                MUT = list(MUT)
                MUT[-(mm_pos+1)] = "T"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_T_REV_MUT"] = "".join(MUT)
            # Replace A (A,G,C)
            if Original == "A":
                temp = list(after)
                # WT REV
                WT = wt + "".join(temp)
                WT = reversecomplement(WT)
                WT = list(WT)
                WT[-(mm_pos+1)] = "A"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_A_REV_WT"] = "".join(WT)
                # MUT REV
                MUT = m + "".join(temp)
                MUT = reversecomplement(MUT)
                MUT = list(MUT)
                MUT[-(mm_pos+1)] = "A"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_A_REV_MUT"] = "".join(MUT)
            if Original == "A":
                temp = list(after)
                # WT REV
                WT = wt + "".join(temp)
                WT = reversecomplement(WT)
                WT = list(WT)
                WT[-(mm_pos+1)] = "G"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_G_REV_WT"] = "".join(WT)
                # MUT REV
                MUT = m + "".join(temp)
                MUT = reversecomplement(MUT)
                MUT = list(MUT)
                MUT[-(mm_pos+1)] = "G"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_G_REV_MUT"] = "".join(MUT)
            if Original == "A":
                temp = list(after)
                # WT REV
                WT = wt + "".join(temp)
                WT = reversecomplement(WT)
                WT = list(WT)
                WT[-(mm_pos+1)] = "C"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_C_REV_WT"] = "".join(WT)
                # MUT REV
                MUT = m + "".join(temp)
                MUT = reversecomplement(MUT)
                MUT = list(MUT)
                MUT[-(mm_pos+1)] = "C"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_C_REV_MUT"] = "".join(MUT)
            # Replace T (T,G,C)
            if Original == "T":
                temp = list(after)
                # WT REV
                WT = wt + "".join(temp)
                WT = reversecomplement(WT)
                WT = list(WT)
                WT[-(mm_pos+1)] = "T"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_T_REV_WT"] = "".join(WT)
                # MUT REV
                MUT = m + "".join(temp)
                MUT = reversecomplement(MUT)
                MUT = list(MUT)
                MUT[-(mm_pos+1)] = "T"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_T_REV_MUT"] = "".join(MUT)
            if Original == "T":
                temp = list(after)
                # WT REV
                WT = wt + "".join(temp)
                WT = reversecomplement(WT)
                WT = list(WT)
                WT[-(mm_pos+1)] = "G"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_G_REV_WT"] = "".join(WT)
                # MUT REV
                MUT = m + "".join(temp)
                MUT = reversecomplement(MUT)
                MUT = list(MUT)
                MUT[-(mm_pos+1)] = "G"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_G_REV_MUT"] = "".join(MUT)
            if Original == "T":
                temp = list(after)
                # WT REV
                WT = wt + "".join(temp)
                WT = reversecomplement(WT)
                WT = list(WT)
                WT[-(mm_pos+1)] = "C"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_C_REV_WT"] = "".join(WT)
                # MUT REV
                MUT = m + "".join(temp)
                MUT = reversecomplement(MUT)
                MUT = list(MUT)
                MUT[-(mm_pos+1)] = "C"
                ART_primers[seq_ID + "_" + str(mm_pos) + "_C_REV_MUT"] = "".join(MUT)
    return ART_primers

ART_primers_REV = REVERSE_TEMPLATES(after, wt, m, position_mismatch, seq_ID)
print("ART_primers_REV:")
print(ART_primers_REV)

##################################################################################################################################

# Just to test the concept

## 1MM (Keep at 60°C)
temps_AS = {}
for key, value in AS_primers.items():
    myseq = Seq(value)
    match = mt.Tm_NN(myseq, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=250, dnac2=250, Na=50, K=0, Tris=75, Mg=3, dNTPs=1.2)
    # get the temperature to be around 60
    while not match <= 60.5:
        myseq = Seq(myseq[1:])
        # stop if the primer gets to short
        if len(myseq) <= 16:
            break
        match = mt.Tm_NN(myseq, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=250, dnac2=250, Na=50, K=0, Tris=75, Mg=3, dNTPs=1.2)
    temps_AS[key] = [str(myseq), '%0.3f' % match, len(myseq)]

print("temps_AS:")   
print(temps_AS)

# Primer 3 should then find a complemantary primer

## FORWARD PRIMERS (Ideally delta of 5°C so around 55°C)
temps_FWD = {}
for key, value in ART_primers_FWD.items():
    myseq = Seq(value)
    original_FWD_WT = Seq(before+wt)
    original_FWD_MUT = Seq(before+m)
    if key.split("_")[-1] == "WT":
        template_2MM = original_FWD_MUT.complement()
        mismatch_TM = mt.Tm_NN(myseq, c_seq=template_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=250, dnac2=250, Na=50, K=0, Tris=75, Mg=3, dNTPs=1.2)
        while not mismatch_TM <= 55.5:
            myseq = Seq(myseq[1:])
            template_2MM = Seq(template_2MM[1:])
            # stop if the primer gets to short
            if len(myseq) <= 16:
                break
            mismatch_TM = mt.Tm_NN(myseq, c_seq=template_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=250, dnac2=250, Na=50, K=0, Tris=75, Mg=3, dNTPs=1.2)
        temps_FWD[key] = [str(myseq), '%0.3f' % mismatch_TM, len(myseq)]
    elif key.split("_")[-1] == "MUT":
        template_2MM = original_FWD_WT.complement()
        mismatch_TM = mt.Tm_NN(myseq, c_seq=template_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=250, dnac2=250, Na=50, K=0, Tris=75, Mg=3, dNTPs=1.2)
        while not mismatch_TM <= 55.5:
            myseq = Seq(myseq[1:])
            template_2MM = Seq(template_2MM[1:])
            # stop if the primer gets to short
            if len(myseq) <= 16:
                break
            mismatch_TM = mt.Tm_NN(myseq, c_seq=template_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=250, dnac2=250, Na=50, K=0, Tris=75, Mg=3, dNTPs=1.2)
        temps_FWD[key] = [str(myseq), '%0.3f' % mismatch_TM, len(myseq)]


print("temps_FWD:")
print(temps_FWD)


## REVERSE PRIMERS (Ideally delta of 5°C so around 55°C)
temps_REV = {}
for key, value in ART_primers_REV.items():
    myseq = Seq(value)
    original_REV_WT = Seq(wt+after).reverse_complement()
    original_REV_MUT = Seq(m+after).reverse_complement()
    if key.split("_")[-1] == "WT":
        template_2MM = original_REV_MUT.complement()
        mismatch_TM = mt.Tm_NN(myseq, c_seq=template_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=250, dnac2=250, Na=50, K=0, Tris=75, Mg=3, dNTPs=1.2)
        while not mismatch_TM <= 55.5:
            myseq = Seq(myseq[1:])
            template_2MM = Seq(template_2MM[1:])
            # stop if the primer gets to short
            if len(myseq) <= 16:
                break
            mismatch_TM = mt.Tm_NN(myseq, c_seq=template_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=250, dnac2=250, Na=50, K=0, Tris=75, Mg=3, dNTPs=1.2)
        temps_REV[key] = [str(myseq), '%0.3f' % mismatch_TM, len(myseq)]
    elif key.split("_")[-1] == "MUT":
        template_2MM = original_REV_WT.complement()
        mismatch_TM = mt.Tm_NN(myseq, c_seq=template_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=250, dnac2=250, Na=50, K=0, Tris=75, Mg=3, dNTPs=1.2)
        while not mismatch_TM <= 55.5:
            myseq = Seq(myseq[1:])
            template_2MM = Seq(template_2MM[1:])
            # stop if the primer gets to short
            if len(myseq) <= 16:
                break
            mismatch_TM = mt.Tm_NN(myseq, c_seq=template_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=250, dnac2=250, Na=50, K=0, Tris=75, Mg=3, dNTPs=1.2)
        temps_REV[key] = [str(myseq), '%0.3f' % mismatch_TM, len(myseq)]

print("temps_REV:")
print(temps_REV)
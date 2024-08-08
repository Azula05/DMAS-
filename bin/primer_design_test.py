#!/usr/bin/env python3
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

sequence = "TGTTATGAACAGCAAATAAAAGAAACTAAAACGATCCTGAGACTTCCACACTGATGCAATCATTCGTCTGTTTCCCATTCTAAACTGTACCCTGTTACTT[A/C]TCCCCTTCCTATGACATGAACTTAACCATAGAAAAGAAGGGGAAAGAAAACATCAAGCGTCCCATAGACTCACCCTGAAGTTCTCAGGATCCACGTGCAG"

before = "TGTTATGAACAGCAAATAAAAGAAACTAAAACGATCCTGAGACTTCCACACTGATGCAATCATTCGTCTGTTTCCCATTCTAAACTGTACCCTGTTACTT"
after = "TCCCCTTCCTATGACATGAACTTAACCATAGAAAAGAAGGGGAAAGAAAACATCAAGCGTCCCATAGACTCACCCTGAAGTTCTCAGGATCCACGTGCAG"
wt = "A"
m = "C"
position_mismatch = "all"
seq_ID = "DMAS-0"

dnac = 250
Na_conc = 50
K_conc = 0
Tris_conc = 75
Mg_conc = 3
dNTPs_conc = 1.2
Single_MM_Tm = 55
Zero_MM_Tm = 60
"""
Creates templates from pre-processed sequence as input for primer3:
    - 2 AS primer (forward, reverse) => 1MM
    - 3 possible artificial mismatch positions (2, 3, 4) => create long primers for temperature and primer validation
"""


def templ_generation_exchange(before, after, wt, m, seq_ID):
    # A) Allele specific primers (1MM)
    # A good length for PCR primers is generally around 18-30 bases

    # Function to get the reverse complement of a sequence
    def reversecomplement(seq):
        my_seq = Seq(seq)
        complement = my_seq.reverse_complement()
        return str(complement)
    
    # A) FORWARD templates
    temp_FWD_WT = before + wt
    temp_FWD_MUT = before + m

    # B) REVERSE templates
    REV_WT = wt + after
    temp_REV_WT = reversecomplement(REV_WT)
    REV_MUT = m + after
    temp_REV_MUT = reversecomplement(REV_MUT)
    templates = {}
    templates[seq_ID + "_FWD_WT"] = temp_FWD_WT
    templates[seq_ID + "_FWD_MUT"] = temp_FWD_MUT
    templates[seq_ID + "_REV_WT"] = temp_REV_WT
    templates[seq_ID + "_REV_MUT"] = temp_REV_MUT
    return(templates)



templates = templ_generation_exchange(before, after, wt, m, seq_ID)


####################################################################################################################################################################

# Match primers needed to start from later
"""
Shorten the templates from 35 to minimum 16 to be as close to the
zero MM tm as possible (match). The remaining sequence is the perfect match primer.
In the library: ID, primer, Tm_0MM, length
"""

def match_primers(myseq, match_TM, dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc):
    myseq = Seq(myseq[-36:])
    match_TM = mt.Tm_NN(myseq, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
    absolute_difference = 2000000
    while not match_TM <= (float(Zero_MM_Tm)-2):
        previous_difference = absolute_difference
        # stop if the primer gets to short
        if len(myseq) < 16:
            break
        else:
            myseq = Seq(myseq[1:])
        match_TM = mt.Tm_NN(myseq, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
        absolute_difference = abs(match_TM - float(Zero_MM_Tm))
        if absolute_difference < previous_difference:
            best_seq = myseq
            best_tm = match_TM
    return str(best_seq),'%0.3f' % best_tm, len(best_seq)

# initiate the dictionary
primers = {}
# WT on WT perfect match
## FORWARD
myseq = templates[seq_ID + "_FWD_WT"]
seq, Tm, length = match_primers(myseq, Zero_MM_Tm, dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc)
GC_content = gc_fraction(seq) * 100
primers[seq_ID + "_0MM_F_WT"] = str(seq), float(Tm), float('%0.3f' % GC_content), int(length)

## REVERSE
myseq = templates[seq_ID + "_REV_WT"]
seq, Tm, length = match_primers(myseq, Zero_MM_Tm, dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc)
GC_content = gc_fraction(seq) * 100
primers[seq_ID + "_0MM_R_WT"] = str(seq), float(Tm), float('%0.3f' % GC_content), int(length)
# MUT on MUT perfect match
## FORWARD
myseq = templates[seq_ID + "_FWD_MUT"]
seq, Tm, length = match_primers(myseq, Zero_MM_Tm, dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc)
GC_content = gc_fraction(seq) * 100
primers[seq_ID + "_0MM_F_MUT"] = str(seq), float(Tm), float('%0.3f' % GC_content), int(length)
## REVERSE
myseq = templates[seq_ID + "_REV_MUT"]
seq, Tm, length = match_primers(myseq, Zero_MM_Tm, dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc)
GC_content = gc_fraction(seq) * 100
primers[seq_ID + "_0MM_R_MUT"] = str(seq), float(Tm), float('%0.3f' % GC_content), int(length)

####################################################################################################################################################################

# Double_MM_exchange(primer_templates, seq_ID):
def exchange (seq, position_mismatch):
    exchanged = {}
    if position_mismatch == "all":
        position_mismatch = [2, 3, 4]
    else:
        position = [position_mismatch]
    seq = list(seq)
    for position in position_mismatch:
        position = position +1
        if seq[-position] == "A":
            seq[-position] = "T"
            exchanged[str(position-1) + "T"] = "".join(seq)
            seq[-position] = "C"
            exchanged[str(position-1) + "C"] = "".join(seq)
            seq[-position] = "G"
            exchanged[str(position-1) + "G"] = "".join(seq)
            seq[-position] = "A"
        elif seq[-position] == "T":
            seq[-position] = "A"
            exchanged[str(position-1) + "A"] = "".join(seq)
            seq[-position] = "C"
            exchanged[str(position-1) + "C"] = "".join(seq)
            seq[-position] = "G"
            exchanged[str(position-1) + "G"] = "".join(seq)
            seq[-position] = "T"
        elif seq[-position] == "C":
            seq[-position] = "A"
            exchanged[str(position-1) + "A"] = "".join(seq)
            seq[-position] = "T"
            exchanged[str(position-1) + "T"] = "".join(seq)
            seq[-position] = "G"
            exchanged[str(position-1) + "G"] = "".join(seq)
            seq[-position] = "C"
        elif seq[-position] == "G":
            seq[-position] = "A"
            exchanged[str(position-1) + "A"] = "".join(seq)
            seq[-position] = "T"
            exchanged[str(position-1) + "T"] = "".join(seq)
            seq[-position] = "C"
            exchanged[str(position-1) + "C"] = "".join(seq)
            seq[-position] = "G"
    return exchanged

def Single_MM_Tm_match(myseq, complement, Single_MM_Tm , dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc):
    mismatch_TM = mt.Tm_NN(myseq, c_seq= complement, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
    absolute_difference = 2000000
    while not mismatch_TM <= ( float(Single_MM_Tm)-2):
        previous_difference = absolute_difference
        # stop if the primer gets to short
        if len(myseq) < 16:
            break
        else:
            myseq = Seq(myseq[1:])
            complement = Seq(complement[1:])
        mismatch_TM = mt.Tm_NN(myseq, c_seq=complement, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
        absolute_difference = abs(mismatch_TM - float(Single_MM_Tm))
        if absolute_difference < previous_difference:
            best_seq = myseq
            best_complement = complement
            best_tm = mismatch_TM
    return str(best_seq), str(best_complement),'%0.3f' % best_tm, len(best_seq)

# WT primer on the MUT 
## FORWARD
forward = primers[seq_ID + "_0MM_F_WT"][0]
# change out the 2,3,4 positions
exchanged = exchange (forward, position_mismatch)
# calculate the Tm for the primers
for mutation ,oligo in exchanged.items():
    # 1MM 
    complement_1MM = str(Seq(oligo).complement()[:-1]) + str(Seq(m).complement())
        # change the length of the oligo to match the target Tm as close as possible
    oligo, complement_1MM ,Single_MM_Tm, length = Single_MM_Tm_match(oligo, complement_1MM, Single_MM_Tm, dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc)
    # 2MM
    complement_2MM = Seq(str(templates[seq_ID + "_FWD_MUT"])[-len(oligo):]).complement()
    # 0MM
    match_Tm = '%0.3f' % mt.Tm_NN(oligo, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
    double_MM_TM = '%0.3f' %  mt.Tm_NN(oligo, c_seq=complement_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
    # GC content
    GC_content = '%0.2f' % (gc_fraction(oligo) * 100)
    # mismatch delta
    mismatch_delta = '%0.3f' % abs(float(Single_MM_Tm) - float(double_MM_TM))
    # add to the dictionary with the primers
    primers[seq_ID + "_" + mutation + "_F_WT"] = oligo, float(match_Tm), float(Single_MM_Tm), float(double_MM_TM), float(mismatch_delta), float(GC_content), len(oligo)

## REVERSE
reverse = primers[seq_ID + "_0MM_R_WT"][0]
# change out the 2,3,4 positions
exchanged = exchange (reverse, position_mismatch)
# calculate the Tm for the primers
for mutation ,oligo in exchanged.items():
    # 1MM 
    complement_1MM = str(Seq(oligo).complement()[:-1]) + str(Seq(m).complement())
        # change the length of the oligo to match the target Tm as close as possible
    oligo, complement_1MM ,Single_MM_Tm, length = Single_MM_Tm_match(oligo, complement_1MM, Single_MM_Tm, dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc)
    # 2MM
    complement_2MM = Seq(str(templates[seq_ID + "_REV_MUT"])[-len(oligo):]).complement()
    # 0MM
    match_Tm = '%0.3f' % mt.Tm_NN(oligo, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
    double_MM_TM = '%0.3f' %  mt.Tm_NN(oligo, c_seq=complement_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
    # GC content
    GC_content = '%0.2f' % (gc_fraction(oligo) * 100)
    # mismatch delta
    mismatch_delta = '%0.3f' % abs(float(Single_MM_Tm) - float(double_MM_TM))
    # add to the dictionary with the primers
    primers[seq_ID + "_" + mutation + "_R_WT"] = oligo, float(match_Tm), float(Single_MM_Tm), float(double_MM_TM), float(mismatch_delta), float(GC_content), len(oligo)

# MUT primer on the WT
# ## FORWARD
forward = primers[seq_ID + "_0MM_F_MUT"][0]
# change out the 2,3,4 positions
exchanged = exchange (forward, position_mismatch)
# calculate the Tm for the primers
for mutation ,oligo in exchanged.items():
    # 1MM 
    complement_1MM = str(Seq(oligo).complement()[:-1]) + str(Seq(m).complement())
        # change the length of the oligo to match the target Tm as close as possible
    oligo, complement_1MM ,Single_MM_Tm, length = Single_MM_Tm_match(oligo, complement_1MM, Single_MM_Tm, dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc)
    # 2MM
    complement_2MM = Seq(str(templates[seq_ID + "_FWD_WT"])[-len(oligo):]).complement()
    # 0MM
    match_Tm = '%0.3f' % mt.Tm_NN(oligo, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
    double_MM_TM = '%0.3f' %  mt.Tm_NN(oligo, c_seq=complement_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
    # GC content
    GC_content = '%0.2f' % (gc_fraction(oligo) * 100)
    # mismatch delta
    mismatch_delta = '%0.3f' % abs(float(Single_MM_Tm) - float(double_MM_TM))
    # add to the dictionary with the primers
    primers[seq_ID + "_" + mutation + "_F_MUT"] = oligo, float(match_Tm), float(Single_MM_Tm), float(double_MM_TM), float(mismatch_delta), float(GC_content), len(oligo)

## REVERSE
reverse = primers[seq_ID + "_0MM_R_MUT"][0]
# change out the 2,3,4 positions
exchanged = exchange (reverse, position_mismatch)
# calculate the Tm for the primers
for mutation ,oligo in exchanged.items():
    # 1MM 
    complement_1MM = str(Seq(oligo).complement()[:-1]) + str(Seq(m).complement())
        # change the length of the oligo to match the target Tm as close as possible
    oligo, complement_1MM ,Single_MM_Tm, length = Single_MM_Tm_match(oligo, complement_1MM, Single_MM_Tm, dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc)
    # 2MM
    complement_2MM = Seq(str(templates[seq_ID + "_REV_WT"])[-len(oligo):]).complement()
    # 0MM
    match_Tm = '%0.3f' % mt.Tm_NN(oligo, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
    double_MM_TM = '%0.3f' %  mt.Tm_NN(oligo, c_seq=complement_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
    # GC content
    GC_content = '%0.2f' % (gc_fraction(oligo) * 100)
    # mismatch delta
    mismatch_delta = '%0.3f' % abs(float(Single_MM_Tm) - float(double_MM_TM))
    # add to the dictionary with the primers
    primers[seq_ID + "_" + mutation + "_R_MUT"] = oligo, float(match_Tm), float(Single_MM_Tm), float(double_MM_TM), float(mismatch_delta), float(GC_content), len(oligo)

# drop the 0MM primers
primers.pop(seq_ID + "_0MM_F_WT")
primers.pop(seq_ID + "_0MM_R_WT")
primers.pop(seq_ID + "_0MM_F_MUT")
primers.pop(seq_ID + "_0MM_R_MUT")

print(primers)
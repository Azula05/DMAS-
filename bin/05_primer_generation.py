#!/usr/bin/env python3
"""
This script generates the specific primers from the given input.
It creates a file with a table containing the following columns:
    - Primer ID
    - Primer sequence
    - Match Tm
    - Single mismatch Tm
    - Double mismatch Tm
    - MM_delta
    - GC content
    - Length
"""
# ARGUMENT HANDLING
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to main snp script')
parser.add_argument("-s",   nargs=1, required=True,     help="bed file with SNPs")
parser.add_argument("-t",   nargs=1, required=True,     help="seq-file with template")
parser.add_argument("-u",   nargs=1, required=True,     help="file with secondary structures")
parser.add_argument("-p",   nargs=1, required=True,     help="Position of the mismatch eighter all, 2, 3 or 4 postions to 5' end")
parser.add_argument("-dnac", nargs=1, required=True,    help="concetration of the oligos µM")
parser.add_argument("-Na",  nargs=1, required=True,     help="Na concentration mM")
parser.add_argument("-K",   nargs=1, required=True,     help="K concentration mM")
parser.add_argument("-Tris", nargs=1, required=True,    help="Tris concentration mM")
parser.add_argument("-Mg",  nargs=1, required=True,     help="Mg concentration mM")
parser.add_argument("-dNTPs", nargs=1, required=True,   help="dNTPs concentration mM")
parser.add_argument("-S",   nargs=1, required=True,     help="Single mismatch optimal Tm")

args = parser.parse_args()

# Input file
seq_file = args.t[0]
# Files with special positions
SNP_file = args.s[0]
sec_str_file = args.u[0]
# Parameters
position_mismatch = args.p[0]
dnac = float(args.dnac[0])
Na_conc = float(args.Na[0])
K_conc = float(args.K[0])
Tris_conc = float(args.Tris[0])
Mg_conc = float(args.Mg[0])
dNTPs_conc = float(args.dNTPs[0])
Single_MM_Tm_goal = float(args.S[0])

############################################################################################################
#################################    FUNCTION: SPLIT SEQUENCE    ###########################################
############################################################################################################
# get all the info from the sequence input file
seq_info = open(args.t[0]).readline().rstrip()
seq, pos, targeted_snp_pos = seq_info.split('\t')
chrom, start, end = pos.replace(':', '-').split('-')
start = int(start)
end = int(end)
targeted_snp_pos = int(targeted_snp_pos)
seq_length = len(seq) - 4
seq_ID = os.path.splitext(args.t[0])[0]
seq_ID = seq_ID.split("/")[-1]

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

# get separate parts of input sequence
before, after, length, wt, m, snp_type = split_up(seq)

############################################################################################################
##########################################    EXCHANGE    ##################################################
############################################################################################################
if snp_type == "exchange":

    # A) GENERATE TEMPLATE SEQUENCES

    ## Function to get the reverse complement of a sequence
    def reversecomplement(seq):
        my_seq = Seq(seq)
        complement = my_seq.reverse_complement()
        return str(complement)
    
    ## FORWARD templates
    temp_FWD_WT = before + wt
    temp_FWD_MUT = before + m

    ## REVERSE templates
    REV_WT = wt + after
    temp_REV_WT = reversecomplement(REV_WT)
    REV_MUT = m + after
    temp_REV_MUT = reversecomplement(REV_MUT)
    ## save templates in dictionary
    templates = {}
    templates[seq_ID + "_FWD_WT"] = temp_FWD_WT
    templates[seq_ID + "_FWD_MUT"] = temp_FWD_MUT
    templates[seq_ID + "_REV_WT"] = temp_REV_WT
    templates[seq_ID + "_REV_MUT"] = temp_REV_MUT

    #=======================================================================================================#

    # B) GENERATE THE SPECIFIC PRIMERS

    primers = {}
    ## Function to exchange nucleotides
    def exchange (seq, position_mismatch):
        exchanged = {}
        ### Which positions need to be exchanged
        if position_mismatch == "all":
            position_mismatch = [2, 3, 4]
        else:
            position = [position_mismatch]
        ### Make a list of the sequence
        seq = list(seq)
        ### Exchange the nucleotides
        for position in position_mismatch:
            position = position +1  # 0-based to 1-based (because 0 is the position of interest)
            ### A -> T, C, G
            if seq[-position] == "A":
                seq[-position] = "T"
                exchanged[str(position-1) + "T"] = "".join(seq)
                seq[-position] = "C"
                exchanged[str(position-1) + "C"] = "".join(seq)
                seq[-position] = "G"
                exchanged[str(position-1) + "G"] = "".join(seq)
                ### Reset the sequence
                seq[-position] = "A"
            ### T -> A, C, G
            elif seq[-position] == "T":
                seq[-position] = "A"
                exchanged[str(position-1) + "A"] = "".join(seq)
                seq[-position] = "C"
                exchanged[str(position-1) + "C"] = "".join(seq)
                seq[-position] = "G"
                exchanged[str(position-1) + "G"] = "".join(seq)
                ### Reset the sequence
                seq[-position] = "T"
            ### C -> A, T, G
            elif seq[-position] == "C":
                seq[-position] = "A"
                exchanged[str(position-1) + "A"] = "".join(seq)
                seq[-position] = "T"
                exchanged[str(position-1) + "T"] = "".join(seq)
                seq[-position] = "G"
                exchanged[str(position-1) + "G"] = "".join(seq)
                ### Reset the sequence
                seq[-position] = "C"
            ### G -> A, T, C
            elif seq[-position] == "G":
                seq[-position] = "A"
                exchanged[str(position-1) + "A"] = "".join(seq)
                seq[-position] = "T"
                exchanged[str(position-1) + "T"] = "".join(seq)
                seq[-position] = "C"
                exchanged[str(position-1) + "C"] = "".join(seq)
                ### Reset the sequence
                seq[-position] = "G"
        return exchanged

    ## Function to shorten the sequence to the optimal Tm for the single mismatch
    def Single_MM_Tm_match(myseq, complement, Single_MM_Tm_goal , dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc):
        mismatch_TM = mt.Tm_NN(myseq, c_seq= complement, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
        absolute_difference = 2000000
        while not mismatch_TM <= (float(Single_MM_Tm_goal)-1.5):
            previous_difference = absolute_difference
            # stop if the primer gets to short
            ##if len(myseq) <= 16:
            ##    break
            #else:
            myseq = Seq(myseq[1:])
            complement = Seq(complement[1:])
            mismatch_TM = mt.Tm_NN(myseq, c_seq=complement, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
            absolute_difference = abs(mismatch_TM - float(Single_MM_Tm_goal))
            if absolute_difference < previous_difference:
                best_seq = myseq
                best_complement = complement
                best_tm = mismatch_TM
        return str(best_seq), str(best_complement),'%0.3f' % best_tm, len(best_seq)
    
    ## WT primer on the MUT 
    ### FORWARD
    forward = templates[seq_ID + "_FWD_WT"][-36:]
    ### change out the 2,3,4 positions
    exchanged = exchange (forward, position_mismatch)
    ### Loop through the exchanged primers
    for mutation ,oligo in exchanged.items():
        #### 1MM
        complement_1MM_W_W = str(Seq(temp_FWD_WT).complement()[-len(oligo):])
            # change the length of the oligo to match the target Tm as close as possible
        oligo, complement_1MM ,Single_MM_Tm, length = Single_MM_Tm_match(oligo, complement_1MM_W_W, Single_MM_Tm_goal, dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc)
        #### 2MM
        complement_2MM = Seq(str(templates[seq_ID + "_FWD_MUT"])[-len(oligo):]).complement()
        #### 0MM
        match_Tm = '%0.3f' % mt.Tm_NN(oligo, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
        double_MM_TM = '%0.3f' %  mt.Tm_NN(oligo, c_seq=complement_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
        #### GC content
        GC_content = '%0.2f' % (gc_fraction(oligo) * 100)
        #### mismatch delta
        mismatch_delta = '%0.3f' % abs(float(Single_MM_Tm) - float(double_MM_TM))
        #### add to the dictionary with the primers
        primers[seq_ID + "_" + mutation + "_F_WT"] = oligo, float(match_Tm), float(Single_MM_Tm), float(double_MM_TM), float(mismatch_delta), float(GC_content), len(oligo)
    #exit("break")
    ### REVERSE
    reverse = templates[seq_ID + "_REV_WT"][-36:]
    # change out the 2,3,4 positions
    exchanged = exchange (reverse, position_mismatch)
    ### Loop through the exchanged primers
    for mutation ,oligo in exchanged.items():
        #### 1MM 
        complement_1MM_W_W = str(Seq(temp_REV_WT).complement())[-len(oligo):]
            # change the length of the oligo to match the target Tm as close as possible
        oligo, complement_1MM ,Single_MM_Tm, length = Single_MM_Tm_match(oligo, complement_1MM_W_W, Single_MM_Tm_goal, dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc)
        #### 2MM
        complement_2MM = Seq(str(templates[seq_ID + "_REV_MUT"])[-len(oligo):]).complement()
        #### 0MM
        match_Tm = '%0.3f' % mt.Tm_NN(oligo, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
        double_MM_TM = '%0.3f' %  mt.Tm_NN(oligo, c_seq=complement_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
        #### GC content
        GC_content = '%0.2f' % (gc_fraction(oligo) * 100)
        #### mismatch delta
        mismatch_delta = '%0.3f' % abs(float(Single_MM_Tm) - float(double_MM_TM))
        #### add to the dictionary with the primers
        primers[seq_ID + "_" + mutation + "_R_WT"] = oligo, float(match_Tm), float(Single_MM_Tm), float(double_MM_TM), float(mismatch_delta), float(GC_content), len(oligo)

    ## MUT primer on the WT
    ### FORWARD
    forward = templates[seq_ID + "_FWD_MUT"][-36:]
    ### change out the 2,3,4 positions
    exchanged = exchange (forward, position_mismatch)
    ### calculate the Tm for the primers
    for mutation ,oligo in exchanged.items():
        #### 1MM 
        complement_1MM_M_M = str(Seq(temp_FWD_MUT).complement())[-len(oligo):]
            #### change the length of the oligo to match the target Tm as close as possible
        oligo, complement_1MM ,Single_MM_Tm, length = Single_MM_Tm_match(oligo, complement_1MM_M_M, Single_MM_Tm_goal, dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc)
        #### 2MM
        complement_2MM = Seq(str(templates[seq_ID + "_FWD_WT"])[-len(oligo):]).complement()
        #### 0MM
        match_Tm = '%0.3f' % mt.Tm_NN(oligo, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
        double_MM_TM = '%0.3f' %  mt.Tm_NN(oligo, c_seq=complement_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
        #### GC content
        GC_content = '%0.2f' % (gc_fraction(oligo) * 100)
        #### mismatch delta
        mismatch_delta = '%0.3f' % abs(float(Single_MM_Tm) - float(double_MM_TM))
        #### add to the dictionary with the primers
        primers[seq_ID + "_" + mutation + "_F_MUT"] = oligo, float(match_Tm), float(Single_MM_Tm), float(double_MM_TM), float(mismatch_delta), float(GC_content), len(oligo)

    ## REVERSE
    reverse = templates[seq_ID + "_REV_MUT"][-36:]
    ### change out the 2,3,4 positions
    exchanged = exchange (reverse, position_mismatch)
    ### calculate the Tm for the primers
    for mutation ,oligo in exchanged.items():
        #### 1MM 
        complement_1MM_M_M = str(Seq(temp_REV_MUT).complement())[-len(oligo):]
            # change the length of the oligo to match the target Tm as close as possible
        oligo, complement_1MM ,Single_MM_Tm, length = Single_MM_Tm_match(oligo, complement_1MM_M_M, Single_MM_Tm_goal, dnac, Na_conc, K_conc,Tris_conc,Mg_conc,dNTPs_conc)
        #### 2MM
        complement_2MM = Seq(str(templates[seq_ID + "_REV_WT"])[-len(oligo):]).complement()
        #### 0MM
        match_Tm = '%0.3f' % mt.Tm_NN(oligo, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
        double_MM_TM = '%0.3f' %  mt.Tm_NN(oligo, c_seq=complement_2MM, nn_table=mt.DNA_NN4, saltcorr=7, dnac1=dnac, dnac2=dnac, Na=Na_conc, K=K_conc, Tris=Tris_conc, Mg=Mg_conc, dNTPs=dNTPs_conc)
        #### GC content
        GC_content = '%0.2f' % (gc_fraction(oligo) * 100)
        #### mismatch delta
        mismatch_delta = '%0.3f' % abs(float(Single_MM_Tm) - float(double_MM_TM))
        #### add to the dictionary with the primers
        primers[seq_ID + "_" + mutation + "_R_MUT"] = oligo, float(match_Tm), float(Single_MM_Tm), float(double_MM_TM), float(mismatch_delta), float(GC_content), len(oligo)

    #=======================================================================================================#

    # C) WRITE PRIMERS TO FILE
    primers_file = open(seq_ID + "_primers.tsv", "a")
    for ID, attributes in primers.items():
        primers_file.write(ID + "\t" + attributes[0] + "\t" + str(attributes[1]) + "\t" + str(attributes[2]) + "\t" + str(attributes[3]) + "\t" + str(attributes[4]) + "\t" + str(attributes[5]) + "\t" + str(attributes[6]) + "\n")
    
############################################################################################################
###########################################    INSERTION   #################################################
############################################################################################################

############################################################################################################
###########################################    DELETION   ##################################################
############################################################################################################

#!/usr/bin/env python3
"""
This script is used to generate files from user input. The user input file is a text file that contains the following
information:
    1. The sequence with the mutation indicated by [].
    2. The chromosome and the coordinates of the sequence.
    Note: make sure the sequence and the coordinates are separated by a tab and not 4 spaces.
    example: 
    TGTTATGAACAGCAAATAAAAGAAACTAAAACGATCCTGAGACTTCCACACTGATGCAATCATTCGTCTGTTTCCCATTCTAAACTGTACCCTGTTACTT[A/C]TCCCCTTCCTATGACATGAACTTAACCATAGAAAAGAAGGGGAAAGAAAACATCAAGCGTCCCATAGACTCACCCTGAAGTTCTCAGGATCCACGTGCAG   chr11:5226403-5226603

If only the sequence is provided, the tool will map the input sequence against a reference genome to determine the information.
If the sequence is provided with the chromosome and coordinates, --coords in the pipeline should be set to 'true'.
"""
import argparse
import os
import re

#######################################################################################################################
# ARGUMENT HANDLING
parser = argparse.ArgumentParser(description="give arguments to handle user input")
parser.add_argument("-c", nargs=1, required=True, help="boolean if sequence coords are manually provided as it should be")
parser.add_argument("-i", nargs=1, required=True, help="user input file")
parser.add_argument("-p", nargs=1, required=True, help="number of cpus to use for process")
parser.add_argument("-b", nargs=1, required=True, help="bowtie index files ; location/basename")

args = parser.parse_args()

coords = args.c[0] # 'true' (coordinates are manually provided) or 'false' (coordinates are not provided)
usr_input_file = open(args.i[0]).readlines()

#######################################################################################################################
####################################### COORDINATES ARE MANUALLY PROVIDED #############################################
#######################################################################################################################

if coords == "true":
    seq_ID_nr = 0
    
    # read input file per 3 lines (Sequence on the first line, coordinates are on every second line, third line is a blank spacer )
    for line in usr_input_file:
        # get the information from the input file
        arguments = line.rstrip().split("\t")
        # check if line is empty
        if len(line.strip()) == 0:
            continue
        # check if the line contains the correct format
        elif len(arguments) != 2:
            line_nr = seq_ID_nr + 1
            message = "Line "+ str(line_nr) + " contains a wrong format. Please check the input file or consider setting coord to false."
            exit(message)
        # extract information from the line
        sequence = arguments[0]
        position = arguments[1]
        # targeted SNP position (0-based)
        seq = list(sequence)
        c = 0
        for nucleotide in seq:
            # the snp comes behind [
            if nucleotide == "[":
                snp_pos = c
            c += 1
        # the separate inpute files
        seq_file = open("DMAS-" + str(seq_ID_nr) + ".txt", "w")
        sequence = sequence.upper()
        seq_file.write(sequence + '\t' + position + '\t' + str(snp_pos) + '\n')
        seq_file.close()
        # Add warning if the sequence is shorter than 150 nt
        if len(sequence) < 150:
            # create a warnings file
            warning_file = open("warning_DMAS-" + seq_ID_nr +".txt", "a")
            warning_file.write("Warning: Input sequence " + str(seq_ID_nr) + " is shorter than 150 nt.\n")
        # increment the sequence ID number
        seq_ID_nr = seq_ID_nr + 1

#######################################################################################################################
##################################### COORDINATES ARE  NOT MANUALLY PROVIDED ##########################################
#######################################################################################################################

elif coords == 'false':

    bowtie_index = args.b[0]
    cpus = args.p[0]

    seq_fasta = open("mapping.fasta","w")

    # loop through input and make a input file for bowtie (fasta) so coordinates can be retrieved
    for i in range(0, len(usr_input_file)):
        # for bowtie2 to recognize the fasta file, > has to precede the sequence
        seq_fasta.write(">DMAS-" + str(i) + "\n")
        seq_fasta.write(usr_input_file[i].rstrip().split("[")[0] + usr_input_file[i].rstrip().split("[")[1].split("/")[0] + usr_input_file[i].rstrip().split("]")[1] + "\n")
    seq_fasta.close()

    #######################################################################################################################
    ########################### use bowtie to map the sequences and get coordinates #######################################
    #######################################################################################################################
    # run bowtie2
    return_code = os.system("bowtie2 -p " + cpus + " --no-hd --no-sq -x " + bowtie_index + " -f " + "mapping.fasta" + " -S mapping.sam")
    # will make debugging easier in nextflow
    if return_code != 0:
        print("Command failed with return code:", return_code)
    # get bowtie results
    sam = open('mapping.sam')
    sam_dict = {}
    fails = []
    
    for line in sam:
        result = line.rstrip().split("\t")
        seq_ID = result[0]
        chrom = result[2]
        start = int(result[3])
        length = result[5]


        # when there is no mapping result
        if length == "*":
            fails.append(seq_ID)

        length = int(re.split(r"(\d+)", length)[1])  # split between numbers and letters
        end = start + length - 1
        sam_dict[seq_ID] = chrom + ":" + str(start) + '-' + str(end)

    # loop through input file again to make a file per input sequence
    usr_input_file = open(args.i[0]).readlines()
    n= 0
    for i in range(0, len(usr_input_file)):
        seq_ID = "DMAS-" + str(i) 
        seq = list(usr_input_file[i])

        # position of the snp
        c = 0
        for i in seq:
            if i == "[":
                snp_pos = c
            c += 1

        seq_file = open(seq_ID + ".txt", "w")
        print(usr_input_file[n])
        seq_file.write(usr_input_file[n].rstrip().upper() + '\t' + sam_dict[seq_ID] + '\t' + str(snp_pos) + '\n')
        seq_file.close()
        n += 1
    if fails != []:
        quit("The following sequences could not be mapped: " + str(fails))
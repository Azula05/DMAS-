#!/usr/bin/env python3
"""
This script will check the generated files and see if the coordinates and the sequence match.
THIS CAN ONLY BE DONE FOR HUMAN GENOME (hg38) AS THE API USED IS SPECIFIC FOR HUMAN GENOME.
If there is a mismatch between the coordinates and the sequence, the mutation will most likely not be in the correct position eighter.
An error message will be given and the script will exit.
"""
# ARGUMENT HANDLING
import argparse
import requests
argparser = argparse.ArgumentParser(description="give arguments to check the coordinates")
argparser.add_argument("-i", nargs=1, required=True, help="generated seq- input file")
args = argparser.parse_args()
file = open(args.i[0]).readlines()

# PARSE THE FILE
for line in file:
    line = line.rstrip().split("\t")
    sequence_input = line[0]
    position = line[1]
    snp_pos = line[2]

    # GET THE SEQUENCE FROM THE COORDINATES
    def fetch_ucsc_sequence(chromosome, start, end, genome='hg38'):
        # Construct the URL for the API request
        url = f"https://api.genome.ucsc.edu/getData/sequence?genome={genome};chrom={chromosome};start={start};end={end}"
        # Make the request to the UCSC API
        response = requests.get(url)
        
        # Raise an exception if the request was unsuccessful
        response.raise_for_status()
        
        # Extract the sequence from the JSON response
        sequence_data = response.json()
        sequence = sequence_data.get("dna", "")
        
        return sequence

    # CHECK IF TEH SEQUENCES MATCH
    ## UCSC
    chromosome = position.split(":")[0]
    start = int(position.split(":")[1].split("-")[0])-1
    end = position.split(":")[1].split("-")[1]
    sequence = fetch_ucsc_sequence(chromosome, start, end)
    sequence = sequence.upper()
    ## INPUT
    sequence_input= sequence_input.split("/")
    left = sequence_input[0].replace("[", "")
    right = sequence_input[1].split("]")[1]
    sequence_input = left + right
    sequence_input = sequence_input.upper()
    ## compare
    if sequence_input != sequence:
        message = "The sequence and the coordinates do not match for "+ str(args.i[0]) +". Please check the input file."
        exit(message)
    ## add message to warnings if the sequences match
    else:
        open("warning.txt", "a").write("The sequence and the coordinates match for "+ str(args.i[0]) +"\n")

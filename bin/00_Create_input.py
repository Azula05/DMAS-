#!/usr/bin/env python3
"""
This script is not included in the pipeline, but added as an extra tool to make the creation of the input file easier.
The script can be ran with the docker file, as the pipeline is ran.
For more information on how to use the script, see the documentation.
Note: this only works for the hg38 genome.
"""

#######################################################################################################################
import argparse
import requests
parser = argparse.ArgumentParser(description="give arguments to handle user input")
parser.add_argument("-i", nargs=1, required=True, help="postitions of interest to create a sequence from")
parser.add_argument("-l", nargs=1, required=False, default=150, help="length of the sequence left and right of the position")

args = parser.parse_args()

# create the input file
input_file = open(args.i[0]).readlines()
try:
    length = int(args.l[0])
except:
    length = 150
output_file = open("Input_file.txt", "w")

#######################################################################################################################
# FUNCTION TO FETCH THE SEQUENCE
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

# get all the information line per line and write the new information
for line in input_file:
    # required information
    line = line.strip().split("\t")
    chromosome = line[0].split(":")[0]
    position = int(line[0].split(":")[1])
    interest = line[1]
    # get the sequence
    start = position - length - 1 # python conversion
    end = position + length
    sequence_left = fetch_ucsc_sequence(chromosome, start, (position-1)) # Do not include the position itself
    sequence_right = fetch_ucsc_sequence(chromosome, position, end) # so there is no extra -1 the postion is not included
    sequence = sequence_left.upper() + interest + sequence_right.upper()
    # write the output
    output_file.write(f"{sequence}\t{chromosome}:{start+1}-{end}\n")

output_file.close()


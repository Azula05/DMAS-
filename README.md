# DMAS
This repository contains the pipeline for designing primers for Double-mismatch Allele-specific qPCR (DMAS-qPCR).

## Installation
### Requirements
The user has the choice of providing the chromosome, start and end position information (--coords). If this information is not provided, the given template is mapped against a reference genome. For this, the user has to [download](https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip) or build the bowtie2 index files of their organism of choice and have them available for the tool (--index_bowtie).

## Usage
The input file should contain the sequence and optionally the sequence coordinates (chr:start-end). It is recommende that the sequence contains at least 150 nucleotides.

The targeted SNP in the sequence should be specified using the following syntax: [N/N]

If coordinates are provided, --coords should be given as an option in the command. Otherwise, the tool will map the input sequence against a reference genome (human by default) to retrieve this information itself. This information is needed to retrieve the common SNPs of the input sequence from a database, because these regions are avoided during primer design.
The coordinates are zero-based, eg. if you give a sequence of 201 nucleotides starting on position 403, you should provide the coordinates 403-603.

Example of input file:
```
TGTTATGAACAGCAAATAAAAGAAACTAAAACGATCCTGAGACTTCCACACTGATGCAATCATTCGTCTGTTTCCCATTCTAAACTGTACCCTGTTACTT[A/C]TCCCCTTCCTATGACATGAACTTAACCATAGAAAAGAAGGGGAAAGAAAACATCAAGCGTCCCATAGACTCACCCTGAAGTTCTCAGGATCCACGTGCAG
chr11:5226403-5226603
```

A help message with available options and arguments is provided by adding the option --help to the command:
```
Example of a typical command for running the pipeline:
$ nextflow run dmas.nf --input template.txt

Required arguments:
    --input                 file with sequence to design primers for
    --outdir		    path to directory where output files should be saved

Optional options:
    --coords		    provide only if input file contains sequence coordinates to skip the mapping step

Optional arguments:
    --settings		    path to file with primer3 settings
    --index_bowtie	    path to bowtie2 index files, including the common file names (excluding .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2)
    --snp_url		    SNP database URL (default: Homo sapiens hg38)
    --spec_filter	    stringency of filtering primer specificity, can be set to 'strict' or 'loose' (default: strict)
    --upfront_filter        When set to 'yes', SNPs and secundary structures are avoided before primer design; when set to 'str', secundary structures are avoided before primer design; when set to 'snp', snp are avoided before primer design; when set to 'no', no filtering before primer design is performed (default: yes)

    Primer3 settings:
    --primer3_diff	    primer min left 3 prime distance (default: 1)
    --primer3_nr	    number of primers to return (default: 20)
    --min_tm		    min melting temp (default: 58)
    --max_tm		    max melting temp (default: 60)
    --opt_tm		    optimal melting temp (default: 59)
    --diff_tm		    melting temp difference between primers (default: 2)
    --min_gc		    min GC content (default: 30)
    --max_gc		    max GC content (default: 80)
    --opt_gc		    optimal GC content (default: 50)
    --amp_min		    min amplicon length (default: 60)
    --amp_max		    max amplicon length (default: 150)
    --mis_lib		    fasta file with mispriming library (default: empty)
    --max_mis_lib	    max allowed weighted similarity with any sequence in mispriming library file (default: 12)
```

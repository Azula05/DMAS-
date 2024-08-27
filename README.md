# DMAS

This repository contains the scripts, docker file and pipeline for designing double mismatch allele specific primers(DMAS) for qPCR (DMAS-qPCR). This tool was developed with the GRCh38 reference genome and human species. However, other species can be used if the correct bowtie2 reference is build. However, the coordinates given in the input will only be checked for the human genome and will only work for hg38. Additionally, the included script 00_Create_input.py will also only work for Hg38; in this case you will have to make the input file yourself.

DMAS primers are created to selectively amplify a specific allele (variant) and obtain an extra mismatch when combining a specific primer with the opposite template (see Figure1). This primer is designed with two intentional mismatches near its 3' end, in addition to the mismatch created by the allele itself. This way, enhance the specificity of the primer is obtained and a more distinct difference can be observed in the qPCR results. These additional mismatches are introduced 2, 3 or 4 positions before (5') of the position of interest or position 0 (see figure 2).

![MM first cycle](https://github.com/user-attachments/assets/5fd1b8e3-1889-4c54-bbfd-5a3e331e5c3f)
Figure 1: primers and template mismatches before PCR.

In Figure 1 all possible combinations are shown before PCR. Both wild type (WT) and mutant (MUT) have a normal primer and a primer where a mismatch is introduced at position -2, -3 or -4. This way a single mismatch (MM) is achieved when combining this special primer with the corresponding template and 2 MM when combined with the other template.

<img src="https://github.com/user-attachments/assets/d61ffd19-47d9-4882-831f-56533c37cfdc" alt="Positions" style="width:300px;"/>
Figure 2: illustration of how the positions are counted.

##  Contents
- [Installation](##Installation)
- [Input](##Input)
	- [Method1](###Method-1:-Location-of-the-position-of-interest)
	- [Method2](###Method-2:-Provide-the-input-directly)
- [Usage](##Usage)
	- [Help](###The-available-parameters)
	- [Master mixes](###Salt-settings)
- [Output](##Output)
- Process (all steps explained)
- Filters: defenities, hoe, elke stap kort uitleggen. wat is SNP etc, structure, specificiteit 
- [Opening a tsv file](##Opening-a-tsv-file-in-excel)

<hr>
## Installation

### Requirements
#### A) Software:
This pipeline can be entirely ran in a docker environment, therefore the only programs that have to be installed locally are: [Docker](https://docs.docker.com/engine/install/) and [Nextflow](https://www.nextflow.io/docs/latest/install.html). 
#### B) Bowtie2 index
This [pipeline](##Explanation-of-the-pipeline) makes use of Bowtie2, to check if the given coordinates match the template and search for off-targets associated with the primer. To be able to do this it requires an index, which can be created by following the steps described below. Bowtie2 does not have to be installed for this process.

- Download the zip file for hg38 in /Assets/GRCh38/Index_bowtie/
```
# Move to the correct file location
cd DMAS/Assets/GRCh38/Index_bowtie/

# Download the files
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip

# Unzip the files
unzip GRCh38_noalt_as.zip
```  
<hr>
## Input
As input this pipeline requires a sequence with a position of interest and the corresponding coordinates in a specific format. Input can be given via 1 of 2 methods. Note that for both methods **the coordinates are zero-based**.

### Method 1: Location of the position of interest
By collecting the chromosomal locations and variant of interest, you can use the included 00_Create_input.py script to make the input for you. This script can be run via the [docker image](https://hub.docker.com/r/oncornalab/dmas) and the output will be called **Input_file.txt**. *Please use this name as input for the pipeline!* By default 150 bases left and right of the template will be picked to create the template. You can change this length by adding the -l option.

To do this your input for this script should look as following (make sure to use a tab!):
	*chrosome:location    [WT/MUT]*
Example:
```
chr11:5226503	[A/C]
chr6:33173384	[A/C]
```

Run the script with the docker image, by replacing [YOUR_INPUT_FILE] with your file ex. variants.txt
```
docker run -v "$PWD":/work -w /work oncornalab/dmas:latest 00_Create_input.py -i /work/[YOUR_INPUT_FILE]
```
or
```
docker run -v "$PWD":/work -w /work oncornalab/dmas:latest 00_Create_input.py -i /work/[YOUR_INPUT_FILE] -l 150
```

### Method 2: Provide the input directly
You can skip the process described in method 1 and create the output from this script yourself if you are working with a sequence which might contain variants such as specific patient sequencing data. To do this your input file should look like the example below, make sure to use tabs instead of four spaces. You can name your file **Input_file.txt**, this way you don't have to change the input file in the nextflow.config file. However, you can give this any name you like as long as it is described in the config file where it is located and how it is called.

Example:
```
TGTTATGAACAGCAAATAAAAGAAACTAAAACGATCCTGAGACTTCCACACTGATGCAATCATTCGTCTGTTTCCCATTCTAAACTGTACCCTGTTACTT[A/C]TCCCCTTCCTATGACATGAACTTAACCATAGAAAAGAAGGGGAAAGAAAACATCAAGCGTCCCATAGACTCACCCTGAAGTTCTCAGGATCCACGTGCAG	chr11:5226403-5226603
TCCTATAACATGGGAGAGTCTTTCTTTTTAGGATTTAGTGGCTCAACAGGTTTGGAAGGACTTGGGAAGGGGAGGCTTCCCTGAGTTCTGTTGCTGGGAGAAGACAGGCAATTCAGTAGGTCAAATGGCAAAGGTAAAACCACAGTGGAAGGCTTCTCTGTGGACAGAGACTGAAGGGTGTGGAGgtat[C/A]GAAGCTGAACGGCAGCAGGGGTGCCTGACAGAATCTCAGCTGCCATCCTCAGGGACTCAGAAGCAGCCTTTTCCGCTTCTGCAGCAATCATCTAGAAAACATGTGACGAAAGCAAAGTGATTGTTCTTCATTCCCTGAAGGCTTCACCACTATGCAGTATGACTCAGCAGACAAGCACTGAGCATCTACTATGTGGCAAGCACGGTTAAGC	chr1:179551171-179551571
```

Make sure your sequence is at least 150 nucleotides long to be able to design primers. The input should be your sequence with [WT/MUT] as your position of interest, followed by a tab and the chromosomal location of this sequence (0-based; UCSC)

<hr>
## Usage

This pipeline is for now only capable of doing **substitutions**, deletions and insertions are not possible in this version of the tool. The targeted positions in the sequence should be specified using the following syntax: [N/N]. With the first base being the wild type and the second one being the variant. Input can eighter be given directly (Method 2) or created from the location of the variant (Method 1). This pipeline is able to handle multiple input lines at once; please not that some processes in this pipeline take time. It may seem like the pipeline is froze but this is not the case. Checking the specificity takes the most time, this can be skipped by turning this option of in the config file. This way the pipeline can be ran faster if you are not interested in off-targets.

While running times might vary depending on hardware these values should give you an indication if you 4 cores available:
- Specificity "loose" or "strict": approximately 5-6 min per input on a personal computer.
- Specificity "off": approximately 15 sec per input when ran on a personal computer.

To run the pipeline you can use the following command (make sure Input_file.txt is in the same folder as DMAS.nf):
```
nextflow run DMAS.nf -profile standard --cpus 3 --input Input_file.txt
```
*This command will instruct nextflow to use the configurations described in the profile standard and that there are 3 threads available to use*. You can change the *standard* to any profile described in the *nextflow.config* file. You can change the number of threads to any number available on your system, however it is better to not use all available threads to keep your system from getting overloaded.

A help message with available options and arguments is provided by adding the option --help to the nextflow run command:
```
nextflow run DMAS.nf --help
```
This will display the help message.

### The available parameters
The easiest way to adapt parameters in the pipeline is to change them in the *nextflow.config* file and just run the pipeline with `nextflow run DMAS.nf -profile standard --cpus 3`. All available arguments are described below.

Required:
- **input**:  path to the file (tab separated) with sequence to design primers for (0-based annotation)
- **index_bowtie**: path to bowtie genome index directory
- **index_bowtie_name**: the basename of the Bowtie index to be searched (the name of any of the index files up to but not including the final .1.ebwt / .rev.1.ebwt / ...)
- **outdir**: path to directory where output files should be saved
- **primer_settings**: path to file with [primer3](https://primer3.org/manual.html) settings (file included in Assets)

Optional pipeline arguments:
- **coords**: true if the input file contains coordinates, false if only sequence is provided (default: true)
- **cpus*:* number of threads to use in bowtie2 (default: 3) => **IMPORTANT: this parameter should be set to the number of CPUs available on the system**
- **snp_url**: SNP database URL (default: http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.bb)

Filters:
- **spec_filter**: stringency of filtering primer specificity, can be set to 'strict', 'loose' or 'off' (default: loose) => **Turn of to run pipeline fast**
- **snp_filter**: when set to 'strict', no common SNPs are allowed in primer sequence; when set to 'loose', no SNPs are allowed on position -2,-3,-4, can be turned off with "off"
- **sec_str_filter**: when set to 'strict', no secondary structure elements are allowed in primer sequence; when set to 'loose', no secondary structure elements are allowed on position -2,-3,-4,-5,-6 can be turned off with "off"
- **validation_filter**: when set to 'strict', both primers need to pass the validation; when set to 'loose', only the specific primer needs to pass the validation; can be turned off with "off"

Primer3 settings:
- **primer3_diff**: the minimum number of base pairs between the 3' ends of any two left primers (see also primer3 PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE)
- **primer3_nr**: number of primers to return (default: 20)
- **min_left_prime** = 1 min left primer distance (default:1)
- **min_tm**: min melting temp (default: 58)
- **max_tm**: max melting temp (default: 62)
- **opt_tm**: optimal melting temp (default: 60)
- **diff_tm**: melting temp difference between primers (default: 2)
- **min_gc**: min GC content (default: 30)
- **max_gc**: max GC content (default: 80)
- **opt_gc**: optimal GC content (default: 50)
- **amp_min**: min amplicon length (default: 50)
- **amp_max**: max amplicon length (default: 150)
- **mis_lib**: fasta file with mispriming library (default: humrep_and_simple; included in Assets)
- **max_mis_lib**: max allowed weighted similarity with any sequence in mispriming library file (default: 12)

Temperature prediction:

- **dnac**: DNA concentration nM (default: 250)
- **na**: Na+ concentration mM (default: 50)
- **k**: K+ concentration mM (default: 0)
- **tris**: Tris concentration mM (default: 75)
- **mg**: Mg2+ concentration mM (default: 3)
- **dNTPs**: dNTPs concentration mM (default: 1.2)
- **position**: position to incorporate the mismatch see Figure 1; options: 2,3,4 or 'all' (default: all)
- **single_MM_Tm**: ideal melting temperature for single mismatches °C (default: 55)

### Salt settings
Above the the default values for the salt settings are displayed under temperature prediction. In the following part different options are given corresponding to different master mixes.

| Name                                                                       | Identifier | Mg++ (mM) | Na+ (mM) | dNTPs (mM) |
| -------------------------------------------------------------------------- | ---------- | --------- | -------- | ---------- |
| Bio-Rad ddPCR Supermix for Probes (No dUTP)                                | #186-3023  | 3.8       | 50       | 0.8        |
| Roche Digital LightCycler 5× DNA Master                                    | /          | 4.5       | 100      | 1.2        |
| Roche Digital LightCycler 5× RNA Master                                    | /          | 3.2       | 100      | 1.5        |
| Sniper 2× dPCR EvaGreen Master Mix（Rox）                                    | #RT002096A | 3.6       | 60       | 0.8        |
| Sniper 2× dPCR Probe Master Mix Plus (Cy5.5)                               | #RT017096A | 3.6       | 60       | 1.0        |
| Sniper 4× Probe Master Mix (Cy5.5)                                         | #RT020096A | 3.8       | 80       | 1.2        |
| Sniper 5× One-step RT-dPCR Probe Super Mix（Cy5.5)                          | #RT018096A | 3.6       | 80       | 1.6        |
| Qiagen QIAcuity Probe PCR Kit*                                             | /          | 6.3       | 13.3     | 1          |
| Qiagen QIAcuity OneStep Advanced Probe Kit*                                | /          | 9         | 50       | 2          |
| Thermofisher: Absolute Q(tm) Universal DNA Digital PCR Master Mix (5X)     | /          | 4.2       | 40       | 1.25       |
| Thermofisher: Absolute Q(tm) Universal DNA Digital PCR Master Mix (5X)<br> | #A52490    | 5.5       | 50       | 1.25       |
| Thermofisher: Absolute Q(tm) Universal DNA Digital PCR Master Mix (5X)<br> | #A55146    | 5.15      | 37.5     | 1.067      |
| Stilla Technologies (all Crystal Digital PCR mastermixes)                  | /          | 5.5       | 60       | 1.4        |

\* subtract 3-4 °C based on additives

## Output
The output in your chosen output directory should look like this:
```
Output/
├── common_alternatives
│   ├── DMAS-0_common_alternatives.txt
│   └── DMAS-1_common_alternatives.txt
├── DMAS-0_primers.tsv
├── DMAS_log.txt
├── Filtered
│   ├── DMAS-0_filtered_loose.tsv
│   ├── DMAS-0_filtered_strict.tsv
│   ├── DMAS-1_filtered_loose.tsv
│   └── DMAS-1_filtered_strict.tsv
└── warnings
    ├── warningDMAS-0.txt
    └── warningDMAS-1.txt
```

- DMAS-0 = the first input line from your original input file. This will increase accordingly to the file. This number is an ID given to the input to seperate them.

- **Common alternatives**: In this folder text files are located which provide an alternative to the common primer used in the main table DMAS-\*\_primers.tsv. In this file all other primers that where found by primer3 are listed these can be used to replace the common primers in the tables if needed. However, it is recommended to check these primers with [primer3plus](https://www.primer3plus.com/) if you do.

- **Filtered**: In this folder 2 files per input are proved, a loose file which gives the filtered output with all filters set to loose and a strict file with all filters set to strict. Since the tool is still under development filters cannot be selected only turned off. The filter parameters are still being evaluated and updated.

- **warnings**: In this folder all warning files are stored. These files should only contain a confirmation that the coordinates match. If not something is wrong with the input.

- **DMAS-\*\_primers.tsv**: This file contains the full table with all unfiltered information per input line. This can be loaded into a program like excel (See: Opening a tsv file in excel) where you can filter and evaluate the data according to your liking. The file contains the following columns:
	- Name
	- Specific_primer
	- Match_Tm
	- Single_MM_Tm
	- Double_MM_Tm
	- MM_delta GC%
	- Lenght
	- Common_primer
	- Match_Tm_common
	- GC%\_common
	- Length_common
	- SNPs_FWD
	- Sec_str_FWD
	- SNPs_REV
	- Sec_str_REV
	- DeltaG_templat
	- FWD_validation
	- REV_validation
	- Amplicon
	- Amp_length
	- Predicted_structure
	- Amplicon_delta_G
	- Forward_specificity
	- Reverse_specificity
	- Specificity_filter_loose
	- Specificity_filter_strict
	- SNP_filter_loose
	- SNP_filter_strict
	- Sec_str_filter_loose
	- Sec_str_filter_strict
	- Validation_filter_loose
	- Validation_filter_strict

- **DMAS_log.txt**: This log file shows how many primer pairs were lost during the different [filter steps](##Filters) and how many remain in the end. The following columns are present:
	- **Name**: ID ex. DMAS-0
	- **Total**: Total amount of sprimer pairs at the start
	- **Specificity_loose**: Primer pairs lost with specificity filter set to loose
	- **Specificity_strict**: Primer pairs lost with specificity filter set to strict
	- **SNP_loose**: Primer pairs lost with SNP filter set to loose
	- **SNP_strict**: Primer pairs lost with SNP filter set to strict
	- **Sec_str_loose**: Primer pairs lost based on template secondary structure with filter set to loose
	- **Sec_str_strict**: Primer pairs lost based on template secondary structure with filter set to strict
	- **Sec_str_amp_loose**: Primer pairs lost based on amplicon secondary structure with filter set to loose
	- **Sec_str_amp_strict**: Primer pairs lost based on amplicon secondary structure with filter set to strict
	- **Validation_loose**: Primer pairs lost based on primer3 primer pair validation with filter set to loose
	- **Validation_strict**: Primer pairs lost based on primer3 primer pair validation with filter set to loose



## Explanation of the pipeline

## Filters
recommended to turn validation off => too stringent evaluate yourself
## Opening a tsv file in excel

1) Open excel
2) Press Open
![Pasted image 20240827154102](https://github.com/user-attachments/assets/44856f00-9276-4408-a4be-403522a92266)
3) Browse to the file location (output folder)
![Pasted image 20240827154147](https://github.com/user-attachments/assets/c7ddf87a-38ee-45f4-83be-2eb8430d07ba)
4) Make sure the dropdown in the bottom right is set to: "All Files (\*.\*)"
![Pasted image 20240827154211](https://github.com/user-attachments/assets/a7faa1d9-9406-4a9d-9a60-eea87f743729)
5) Open the file
6) A pop-up will appear, select delimited and next
![Pasted image 20240827154255](https://github.com/user-attachments/assets/28e33377-7908-4522-9eb2-6647462be485)
7) Select Tab
![Pasted image 20240827154310](https://github.com/user-attachments/assets/4e9bf009-214c-421d-bae4-35b6f5984338)
8) Next
9) Select general
10) Press Finish
![Pasted image 20240827154344](https://github.com/user-attachments/assets/2548a61a-a8eb-40e8-8b80-f0d9d7546fa2)

# DMAS

This repository contains the scripts, docker file and pipeline for designing double mismatch allele specific primers(DMAS) for qPCR (DMAS-qPCR). This tool was developed with the GRCh38 reference genome and human species. However, other species can be used if the correct bowtie2 reference is build. However, the coordinates given in the input will only be checked for the human genome and will only work for hg38. Additionally, the included script 00_Create_input.py will also only work for Hg38; in this case you will have to make the input file yourself.

DMAS primers are created to selectively amplify a specific allele (variant) and obtain an extra mismatch when combining a specific primer with the opposite template (see Figure1). This primer is designed with two intentional mismatches near its 3' end, in addition to the mismatch created by the allele itself. This way, enhance the specificity of the primer is obtained and a more distinct difference can be observed in the qPCR results. These additional mismatches are introduced 2, 3 or 4 positions before (5') of the position of interest or position 0 (see figure 2).

<img src="https://github.com/user-attachments/assets/5fd1b8e3-1889-4c54-bbfd-5a3e331e5c3f" alt="MM explanation" style="width:500px;"/>

#Figure1: primers and template mismatches before PCR.

In Figure 1 all possible combinations are shown before PCR. Both wild type (WT) and mutant (MUT) have a normal primer and a primer where a mismatch is introduced at position -2, -3 or -4. This way a single mismatch (MM) is achieved when combining this special primer with the corresponding template and 2 MM when combined with the other template.

<img src="https://github.com/user-attachments/assets/d61ffd19-47d9-4882-831f-56533c37cfdc" alt="Positions" style="width:500px;"/>

#Figure2: illustration of how the positions are counted.

Throughout the documentation an example will be shown in each steps if you want to familiarize yourself with the pipeline, you can try this example first before running your own sequences. These examples will be marked with <span style="color:green"> green text</span>.


##  Contents
- [Installation](#Installation)
- [Input](#Input)
	- [Method1](#Method-1-Location-of-the-position-of-interest)
	- [Method2](#Method-2-Provide-the-input-directly)
- [Usage](#Usage)
	- [Help](#The-available-parameters)
	- [Master mixes](#Salt-settings)
- [Output](#Output)
- Process (all steps explained)
- Filters: defenities, hoe, elke stap kort uitleggen. wat is SNP etc, structure, specificiteit 
- [Opening a tsv file](#Opening-a-tsv-file-in-excel)

<hr>


## Installation

### Requirements
#### A) Software:
This pipeline can be entirely ran in a docker environment, therefore the only programs that have to be installed locally are: [Docker](https://docs.docker.com/engine/install/) and [Nextflow](https://www.nextflow.io/docs/latest/install.html). 
#### B) Bowtie2 index
This [pipeline](#Explanation-of-the-pipeline) makes use of Bowtie2, to check if the given coordinates match the template and search for off-targets associated with the primer. To be able to do this it requires an index, which can be created by following the steps described below. Bowtie2 does not have to be installed for this process.

- Download the zip file for hg38 in /Assets/GRCh38/Index_bowtie/
```
# Move to the correct file location (in DMAS folder)
cd ./Assets/GRCh38/Index_bowtie/

# Download the files
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip

# Unzip the files
unzip GRCh38_noalt_as.zip
```

<hr>


## Input
As input this pipeline requires a sequence with a position of interest and the corresponding coordinates in a specific format. Input can be given via 1 of 2 methods. Note that for both methods **the coordinates are zero-based**.

### Method 1 Location of the position of interest

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

#### <span style="color:green">Illustration of method 1</span>
Example.txt (in example folder):
```[example.txt]
TGTTATGAACAGCAAATAAAAGAAACTAAAACGATCCTGAGACTTCCACACTGATGCAATCATTCGTCTGTTTCCCATTCTAAACTGTACCCTGTTACTT[A/C]TCCCCTTCCTATGACATGAACTTAACCATAGAAAAGAAGGGGAAAGAAAACATCAAGCGTCCCATAGACTCACCCTGAAGTTCTCAGGATCCACGTGCAG	chr11:5226403-5226603
```

Run the script on the file and set the length of the template to 160 bases on each side:
```
docker run -v "$PWD":/work -w /work oncornalab/dmas:latest 00_Create_input.py -i /work/Example/Example.txt -l 160
```

Expected output for Input_file.txt
```
ATAGTAAAAATTGCGGAGAAGAAAAAAAAAGAAAGCAAGAATTAAACAAAAGAAAACAATTGTTATGAACAGCAAATAAAAGAAACTAAAACGATCCTGAGACTTCCACACTGATGCAATCATTCGTCTGTTTCCCATTCTAAACTGTACCCTGTTACTT[A/C]TCCCCTTCCTATGACATGAACTTAACCATAGAAAAGAAGGGGAAAGAAAACATCAAGCGTCCCATAGACTCACCCTGAAGTTCTCAGGATCCACGTGCAGCTTGTCACAGTGCAGCTCACTCAGTGTGGCAAAGGTGCCCTTGAGGTTGTCCAGGTGAGC	chr11:5226343-5226663
```

### Method 2 Provide the input directly

You can skip the process described in method 1 and create the output from this script yourself if you are working with a sequence which might contain variants such as specific patient sequencing data. To do this your input file should look like the example below, make sure to use tabs instead of four spaces. You can name your file **Input_file.txt**, this way you don't have to change the input file in the nextflow.config file. However, you can give this any name you like as long as it is described in the config file where it is located and how it is called.

Example:
sequence[N/N]sequence	chr:start-end
```
TGTTATGAACAGCAAATAAAAGAAACTAAAACGATCCTGAGACTTCCACACTGATGCAATCATTCGTCTGTTTCCCATTCTAAACTGTACCCTGTTACTT[A/C]TCCCCTTCCTATGACATGAACTTAACCATAGAAAAGAAGGGGAAAGAAAACATCAAGCGTCCCATAGACTCACCCTGAAGTTCTCAGGATCCACGTGCAG	chr11:5226403-5226603
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
*This command will instruct nextflow to use the configurations described in the profile standard and that there are 3 threads available to use*. You can change the *standard* to any profile described in the *nextflow.config* file. You can change the number of threads to any number available on your system, however it is better to not use all available threads to keep your system from getting overloaded. <span style="color:red">Important: do not run the pipeline again with -resume</span>. Since many different processes require output to be calculated NextFlow will rerun preceding processes as well, which will cause certain columns to be overwritten in the table with the wrong input. If you want to run the pipeline again you will have to run the entire pipeline. Since the specificity filter is at the end small mistakes before this will not result in long waiting times.

A help message with available options and arguments is provided by adding the option --help to the nextflow run command:
```
nextflow run DMAS.nf --help
```
This will display the help message.

#### <span style="color:green">Run the pipeline on the example</span>
```
nextflow run DMAS.nf -profile standard --cpus 3 --input Input_file.txt --outdir ./Example/
```

See folder Example for the expected output.


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


<hr>


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

<div></div>


- **Common alternatives**: In this folder text files are located which provide an alternative to the common primer used in the main table DMAS-\*\_primers.tsv. In this file all other primers that where found by primer3 are listed these can be used to replace the common primers in the tables if needed. However, it is recommended to check these primers with [primer3plus](https://www.primer3plus.com/) if you do.

<div></div>


- **Filtered**: In this folder 2 files per input are proved, a loose file which gives the filtered output with all filters set to loose and a strict file with all filters set to strict. Since the tool is still under development filters cannot be selected only turned off. The filter parameters are still being evaluated and updated.

<div></div>


- **warnings**: In this folder all warning files are stored. These files should only contain a confirmation that the coordinates match. If not something is wrong with the input.

<div></div>


- **DMAS-\*\_primers.tsv**: This file contains the full table with all unfiltered information per input line. This can be loaded into a program like excel (See: Opening a tsv file in excel) where you can filter and evaluate the data according to your liking. The file contains the following columns:
	- **Name**: <span style="color:green">DMAS-0_2A_F_WT</span>
		- DMAS-0 = the ID dependent on the input line, where 0 is the first input
		- 2A: at position -2 the original base was replaced with an A
		- F: this means this primer is the forward primer where R would be the reverse primer
		- WT: This means it a wild type primer (MUT for mutant), this means position 0 will be the wild type variant.
	- **Specific_primer**: Sequence of the altered primer depending on F or R this might be the forward or reverse primer. Look at the name to see the alteration.
	- **Match_Tm**: This is the meting temperature predicted for a perfect match. This means that the altered primer matches its template perfectly, this will be the case after PCR. For example altered WT primer on WT template after 1 PCR cycle.
	- **Single_MM_Tm**: This will be the predicted meting temperature for a single mismatch on position -2, -3 or -4. This means that when the altered primer is matched with its corresponding template WT>WT and MUT>MUT only a single mismatch will be present (the altered position) before PCR. This temperature will be as close as possible to the  *single_MM_Tm* parameter.
	- **Double_MM_Tm**: This is the predicted melting temperature when there are 2 mismatches this would be 1 of the 3 possible positions and position 0 (see [Figure 2](#Figure2)). This only happens when WT>MUT or MUT>WT (see [Figure 1](#Figure1))
	- **MM_delta**: Temperature difference between single MM and double MM.
	- **GC%**: GC content of the altered primer.
	- **Lenght**: length of the altered primer.
	- **Common_primer**: Common primer, this primer will both match the WT and MUT template perfectly without mismatches. If F is in the name than this will be the reverse primer, if R is in the name than this will be the forward primer.
	- **Match_Tm_common**: Predicted melting temperature of the common primer which should be as close as possible to the *optimal Tm* paramter you gave in the parameters.
	- **GC%\_common**: GC content of the common primer.
	- **Length_common**: Length of the common primer.
	- **SNPs_FWD**: Single nucleotide polymorphisms found in the forward primer based on the SNP_URL (default: UCSC common SNP database). Example: <span style="color:green">{-7: ('rs7946748', 'G/A,')}</span>
		- -7 position counted back from 0, so 7 position 5' of the position of interest
		- rs identifier
		- The variant [WT/MUT]
	- **Sec_str_FWD**: Secondary structural element predicted in the template that overlap the forward primer. If such structures are found it would mean that the template is less accessible there for the primers. Example: <span style="color:green">[-3, -2, -1]</span>.
		- The first 3 position 5' of the position of interest are part of a secondary structural element in the template. (Note: look at delta_G)
	- **SNPs_REV**: Single nucleotide polymorphisms found in the reverse primer based on the SNP_URL
	- **Sec_str_REV**: Secondary structural element predicted in the template that overlap the reverse primer.
	- **DeltaG_template**: Delta G value associated with the secondary structure prediction of the template. Very low numbers are more likely to occur that higher numbers (close to 0 or positive).
	- **FWD_validation**: Validation of the forward primer by Primer3. Ideally this value should be : "considered 1, ok 1"; however if the primer does not meet the primer3 settings criteria (see above). Then it will say something like GC content failed. Note: if the message displays "considered 1, ok 0" it might be that the primer can still work. Consider checking with another common primer or taking all other information into account before discarding primers with this value. However, it is an indication that this primer might be less ideal. Primer3 failed often means the amplicon is too small.
	- **REV_validation**: Validation of the reverse primer by Primer3 (see FWD_validation).
	- **Amplicon**: Amplicon produced by the primers (including the primers).
	- **Amp_length**: Lenght of this amplicon.
	- **Predicted_structure**: Predicted secondary structure of this amplicon where ( or ) means that that position is involved in a secondary structure element. Note: consider delta G value to see if this likely to happen. This prediction was made with RNAfold (DNA specific settings).
	- **Amplicon_delta_G**: Delta G value of the secondary structure predicted for the amplicon. Lower values are more likely to occur. Close 0 is not likely to occur.
	- **Forward_specificity**: Specificity of the forward primer. Ideally this should be a chromosomal location withing the template region, in this case mismatches can be ignored be cause the primer is on target and the mismatches are intended. If this would not be the case an off-target location is found, which would mean the primer is not specific enough (can only be ignored if sufficient MM are present; see filters).
	- **Reverse_specificity**: Specificity of the reverse primer.
	- **Specificity_filter_loose**: Conclusion tag for this filter would eighter be "NA" of "off" when switched off or "Pass"/"Fail" if on. These tags can be used to filter yourself by loading the file into a program like excel.
	- **Specificity_filter_strict**: Conclusion tag for this filter would eighter be "NA" of "off" when switched off or "Pass"/"Fail" if on. These tags can be used to filter yourself by loading the file into a program like excel.
	- **SNP_filter_loose**: Conclusion tag for this filter would eighter be "NA" of "off" when switched off or "Pass"/"Fail" if on. These tags can be used to filter yourself by loading the file into a program like excel.
	- **SNP_filter_strict**: Conclusion tag for this filter would eighter be "NA" of "off" when switched off or "Pass"/"Fail" if on. These tags can be used to filter yourself by loading the file into a program like excel.
	- **Sec_str_filter_loose**: Conclusion tag for this filter would eighter be "NA" of "off" when switched off or "Pass"/"Fail" if on. These tags can be used to filter yourself by loading the file into a program like excel. 
	- **Sec_str_filter_strict**: Conclusion tag for this filter would eighter be "NA" of "off" when switched off or "Pass"/"Fail" if on. These tags can be used to filter yourself by loading the file into a program like excel.
	- **Validation_filter_loose**: Conclusion tag for this filter would eighter be "NA" of "off" when switched off or "Pass"/"Fail" if on. These tags can be used to filter yourself by loading the file into a program like excel. (This filter is very stringent, consider this one to narrow down the option instead of dropping options)
	- **Validation_filter_strict**: Conclusion tag for this filter would eighter be "NA" of "off" when switched off or "Pass"/"Fail" if on. These tags can be used to filter yourself by loading the file into a program like excel. (This filter is very stringent, consider this one to narrow down the option instead of dropping options)

<div></div>


- **DMAS_log.txt**: This log file shows how many primer pairs were lost during the different [filter steps](#Filters) and how many remain in the end. The following columns are present:
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

<div></div>


### Picking primers
Important when picking primer for your PCR run, <span style="color:red">do not pick complementary primers</span>. If you would choose for example \_2A_F_WT and  \_2T_F_MUT position -2 would contain eighter an A or a T which is complementary and the mismatch WT>MUT or MUT>WT would be lost. Consider options like \_2A_F_WT with \_2G_F_WT or a different position like \_2A_F_WT with \_4G_F_WT. 

<hr>

<div></div>


## Explanation of the pipeline

<div></div>


### PROCESS 1 - splitting input file

In this process the input file is split up into separate files. If the *coord* parameter is set to true the coordinates are provided (should be the case if Method 1 was used for input), ideally this should also be the case if method 2 was used and you provided the sequences yourself. Every line of the input file is considered an input and is split up. If the line does not contain 2 values (sequence and coordinates) the program will notify the user that the coordinates are missing and that the input file should be checked or *coord* should be set to false (possible but not advised). 
From this input a new file will be written named according to the DMAS-id which will be used to keep files from a certain input together during the whole pipeline so they do not mix. This new file includes the sequence, the coordinates and the position of the position of interest within the provided sequence.
If the coordinates are not provided and *coord* is set to false, the pipeline will try to recover the coordinates by mapping the wild type template with [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) against the provided Bowtie2 index. If it fails the program will notify the user which sequence could not be mapped. Note that while this option could be nice to reduce manual input it lacks user control. If an error would get into the sequence or coordinates at this point it will be taken through the whole pipeline.
Lastly, the log file is created and the header line is written to it.

<div></div>


### PROCESS 2 - check coordinates

This process includes a script which will match the provided sequence against the provided coordinates. If the wild type sequence does not match the hg38 sequence retrieved from [api.genome.ucsc.edu](http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=chr1:179551171,179551571) an error will be displayed showing which sequence did not match its coordinates. Otherwise the warning file should only contain a confirmation that the sequence does match. If an error is displayed and you are certain that these coordinates are correct you can make the pipeline ignore this check by changing species to any other value than human. Note that this step is only performed when the species is human therefore only use this option if you are working with the human hg38 genome. This tool was developed with the human genome in mind, however if the bowtie index is available for other species, the pipeline can be ran for other species as well and this step will not be performed.

<div></div>


### PROCESS 3 - Secondary structure of the template

This step uses [RNAfold](https://www.tbi.univie.ac.at/RNA/) by ViennaRNA to predict the minimal free energy structure of the template that was provided to evaluate if they are accessible for the primers. Despite the RNA being in the name, this tool is also capable of working with DNA when configured correctly. The structure is both predicted for the wild type template and the mutant template. The prediction is based on working temperatures of *opt_tm* for DNA and salt correction based around *na*. The output will be written to files containing the dot bracket notation and the delta G value of the prediction. These will be used in later steps when the primers are designed.

<div></div>


### PROCESS 4 - Find common single nucleotide polymorphisms

In the next step, the coordinates are cross-referenced with the common single nucleotide polymorphism database (UCSC). If any SNPs are found these are written to a bed file which will later be consulted to see if they are in the primer sequences are not.

<div></div>



### PROCESS 5 - Primer generation

This part will generate the 36 adapted primers (WT and MUT). First it will check which kind of mutation is at the position of interest. For now, the pipeline can only handle exchange. Then the wild type and mutant templates are generated both the sense and antisense (reverse complement of the sense template) templates are made.
#### Exchange
Next the altered primers (containing a mismatch at position -2, -3 or -4) are made. The script will check which base was originally present at that location in the template and change the position to the other possibilities. Example: if A was present at position -2, this can be changed to T (T>T), C (C>T) and G (G>T). This is done both for the forward and reverse primers of WT and MUT where the 3' end will always be locked on position 0 (the position of interest). Once all possibilities are created the templates which are now 36 bases long will be evaluated for their melting temperature.
For this predication [BioPython](https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html) is used with salt correction based on Owczarzy et al. (2008) and thermodynamic NN values for DNA from SantaLucia & Hicks (2004). Temperature prediction uses concentrations for Mg++, Na+, DNA, K+, Tris, dNTPs.
The sequence is shortened from 36 until the temperature is as close as possible to the *single_MM_Tm* this means that when the adapted sequence WTprimer>WTtemplate is matched only a single mismatch is present and the prediction of its melting temperature is as close to the provided value as possible (default: 55). Once the length of the primer is found the other temperatures: perfect match (WTprimer>WTtemplate and MUTprimer>MUTtemplate) and double mismatch tm (WTprimer>MUTtemplate and MUTprimer>WTtemplate) are predicted as well. Additionally the GC content is calculated and the lengt is counted. Once this is performed for all primers the name, primer, Tms, GC% and length are written to a table.

#### Insertion/deletion
Not possible yet.

<div></div>


### PROCESS 6 - Common primer generation

Since nextflow cannot use input as output the individual input files (example: DMAS-0.tsv) should be converted to a tuple (DMAS-0, DMAS-0.tsv) this way NextFlow can combine all files starting with the correct ID without matching wrong files.
In the next step the most important [Primer3 setting values](#The-available-parameters)  of the Primer3_settings file can be overwritten with user input without altering the file.
This script will first create 2 primer input files: one where the forward primer is already found and we want possible reverse primers and one where the opposite is true. Since the created primers are only used to indicate where the forward or reverse primer is located the first F_WT and first R_WT sequence is used as placeholder. (Mismatches do not mater in this step since the common primer will not have any and will be based on the wild type template).
All provided settings will be written to the input file these include: SEQUENCE_ID, SEQUENCE_TEMPLATE,  FORWARD or REVERSE primer that was already obtained, PRIMER_NUM_RETURN, PRIMER_MIN_THREE_PRIME_DISTANCE, PRIMER_PRODUCT_SIZE_RANGE, 
PRIMER_MIN_TM, PRIMER_MAX_TM, PRIMER_OPT_TM, PRIMER_PAIR_MAX_DIFF_TM, PRIMER_MIN_GC, PRIMER_MAX_GC, PRIMER_OPT_GC_PERCENT, PRIMER_OPT_GC_PERCENT, PRIMER_MISPRIMING_LIBRARY, PRIMER_INTERNAL_MISHYB_LIBRARY, PRIMER_MAX_LIBRARY_MISPRIMING, PRIMER_DNA_CONC, PRIMER_SALT_MONOVALENT, PRIMER_SALT_DIVALENT, PRIMER_DNTP_CONC.
Afterwards, the file containing all SNPs found in the template is considered and the positions containing such an SNP are added to the list of positions to exclude (if SNP_filter is "off", this will not happen). The position of interest is also added as a position to exclude since we do not want this to be in the common primers. Lastly the region to search for primers is limited to eighter the remaining sequence to the right (3') of the altered primer for the common reverse or the remaining sequence to the left (5') of the altered primer for the common forward.
These input files are then given to Primer3 to find common primers, all possible primers are written to files named Common_REV.txt or Common_FWD.txt to be used in later steps.
### PROCESS 7 - Validate primers

In this process the common primers are first added to the table, in general the first primer of the list with found primers is used since this primer has the least penalties. The other common primers are kept in the folder Alterntive_common_primers, if the amplicon is not satisfactory or you would like to change the common primer this is possible by picking one from this list.
Following, validation of the primer pairs is performed and added to the table in new columns. The altered position of the altered primer is considered since Primer3 requires an exact match. Therefore the templates are temporarily changed to match this alteration creating a perfect match. Next validation of the primer pairs is checked with Primer3 check primers (can performed [here](https://www.primer3plus.com/index.html) by setting task to check primers; for more information; make sure to alter the template accordingly). If for any reason Primer3 fails to check the primer pair the second best common primer is considered instead. If this one also fails, "primer3 failed to validate primer pair" will be outputted to the table. Both validation of the forward and reverse primer will be separate values in the output table.
Next (if SNP filter is on) the found SNPs will be checked with the primer locations if SNPs are found within the primers the relative location, rs identifier and SNP will be written to the table. Additionally, should the position of interest be recognised as a common SNP this will be ignored.
Following, the secondary structure of the template will be checked for accessibility (if on). If secondary structural elements are found within the primer locations, their relative location will be written to the table. Following, the amplicon will be created, if this process fails for any reason NA will be written instead. Next, the length of the amplicon (including primers) is calculated and written to the table.
Lastly, all information is written to the table the table is given a new column order. Additionally, the alternative primer lists are stored in a folder *common_alternatives* within the output folder.
### PROCESS 8 - Get secondary structure of the amplicon

This step uses [RNAfold](https://www.tbi.univie.ac.at/RNA/) by ViennaRNA to predict the minimal free energy structure of the amplicon that was created with the primer pair. Despite the RNA being in the name, this tool is also capable of working with DNA when configured correctly. The prediction is based on working temperatures of *opt_tm* for DNA and salt correction based around *na*. The output will be written to the table in dot bracket notation and the delta G value of the prediction.
### PROCESS 9 - Check the specificity of the primers

This script will use Bowtie2 to map the primers against the full genome. If highly similar sequences are found bowtie2 will output these together with the amount of mismatches. The command used in bowtie2 is: bowtie2 --no-hd --xeq --no-sq -X 1000 -N 1 --mp 1,1 --quiet -x  bowtie2_index  -X1000 --very-sensitive -f - --threads. This will be performed for each primer pair and the positions where the primers could match are outputted to the table together with information about the mismatches. Example: chr1:114430787:18=1X2=:1; this string indicates that a match was found on chromosome 1 with starting position 114430787, than there 18 bases with an exact match a mismatch and 2 more bases with an exact match, resulting in a single mismatch (number at the end). Which in this case was an on target match of a primer with an alteration at position -2, which gave the single mismatch.
This process takes a long time because the whole genome is checked with 72 sequences to match. If you want to run the pipeline faster and specificity is not your main concern, you can run the pipeline with the specificity filter of which significantly speeds up the process. However note that whole pipeline needs to be ran again if you want the septicity. 
### PROCESS 10 - Filter the table based on the different filter parameters

The last step of the process creates an output directory for the filtered tables, and copies the full table with all information to the output folder. A file will be created where all filters are set to loose and another where all filters are set to strict. During this process all filters described in the next part are performed an the amount of primer pairs lost during each step are counted and stored within a log file as well as the amount that passed all filters. Evaluate for yourself which filters you find important for your purpose and turn these on by using "loose" or "strict". All filters can be turned of as well by using "off". The script will check if the output from Bowtie2 is outside of the template coordinates and if sufficient mismatches (see below) are present. If there are not sufficient mismatches the primer is not specific enough and will receive a "Fail" tag in the corresponding columns. Next, the SNPs found in the primers will be evaluated according to the filters described below and will receive "Pass" or "Fail" in their corresponding columns. Following, the template secondary structure and amplicon secondary structures are evaluated as well and a "Pass" or "Fail" tag will be written to the corresponding columns. Lastly, the validation results from Primer3 are evaluated and a "Pass" or "Fail" tag will be written to the corresponding columns.
Turning any of the filters off will result in "off" instead of "PASS" or "FAIL", these are regarded as a "PASS" tag and counted as such in the log file. While the secondary structure filter considers both template and amplicon these will not be split up in the table as is the case in the log file. PASS means both structures passed the filter, while FAIL will here always indicate which one failed; the template or the amplicon.

At the end of this process NextFlow will display some general information about the pipeline, like the run time, profile, input and output as well as the docker image used if any.

<hr>


## Filters

#### Specificity
First the output form Bowtie2 is retrieved from the table, if the chromosome matches the one from the input then the start position will be compared to the coordinates of the template. If the start position in within these coordinates the primer automatically passes. Only if there is another match found where this is not the case and the criteria are not met will the primer fail. The following criteria are used for the specificity filter are:

Strict:
- no mismatches in either primer = off-target
- primers with at least 4 mismatches for a single primer = no off-target
- primers with a total of at least 5 mismatches between both primers = no off-target

<div></div>

 Loose:
 - primers with at least 3 mismatches for a single primer = no off-target
 - primers with a total of at least 4 mismatches between both primers = no off-target

<div></div>

| MM primer 1 | MM primer 2 | Sum MM | Loose                                          | Strict                                         |
| ----------- | ----------- | ------ | ---------------------------------------------- | ---------------------------------------------- |
| 0           | 0           | 0      | <span style="color:red">off-target</span>      | <span style="color:red">off-target</span>      |
| 1           | 0           | 1      | <span style="color:red">off-target</span>      | <span style="color:red">off-target</span>      |
| 0           | 1           | 1      | <span style="color:red">off-target</span>      | <span style="color:red">off-target</span>      |
| 2           | 0           | 2      | <span style="color:red">off-target</span>      | <span style="color:red">off-target</span>      |
| 0           | 2           | 2      | <span style="color:red">off-target</span>      | <span style="color:red">off-target</span>      |
| 1           | 1           | 2      | <span style="color:red">off-target</span>      | <span style="color:red">off-target</span>      |
| 2           | 1           | 3      | <span style="color:red">off-target</span>      | <span style="color:red">off-target</span>      |
| 1           | 2           | 3      | <span style="color:red">off-target</span>      | <span style="color:red">off-target</span>      |
| 3           | 0           | 3      | <span style="color:green">no off-target</span> | <span style="color:red">off-target</span>      |
| 0           | 3           | 3      | <span style="color:green">no off-target</span> | <span style="color:red">off-target</span>      |
| 1           | 3           | 4      | <span style="color:green">no off-target</span> | <span style="color:red">off-target</span>      |
| 3           | 1           | 4      | <span style="color:green">no off-target</span> | <span style="color:red">off-target</span>      |
| 2           | 2           | 4      | <span style="color:green">no off-target</span> | <span style="color:red">off-target</span>      |
| 4           | 0           | 4      | <span style="color:green">no off-target</span> | <span style="color:green">no off-target</span> |
| 0           | 4           | 0      | <span style="color:green">no off-target</span> | <span style="color:green">no off-target</span> |
| 2           | 3           | 5      | <span style="color:green">no off-target</span> | <span style="color:green">no off-target</span> |
| 3           | 2           | 5      | <span style="color:green">no off-target</span> | <span style="color:green">no off-target</span> |
| 4           | 1           | 5      | <span style="color:green">no off-target</span> | <span style="color:green">no off-target</span> |
| 1           | 4           | 5      | <span style="color:green">no off-target</span> | <span style="color:green">no off-target</span> |
| 3           | 3           | 6      | <span style="color:green">no off-target</span> | <span style="color:green">no off-target</span> |

<div></div>
#### SNP filter
If 0 SNPs are found in the primer, it automatically passes the filter. If there are SNPs found the relative position within the primer is considered. This is only considered for the altered primer, if common SNPs are present in the common primer this will be reported but not filtered on The criteria used in this filter are:

Loose:
- If the SNPs are not on position -2, -3 or -4 (they can never be on position 0). Than the primers fail.
- If not the primer passes
Stict:
- The primer pair only passes when no SNPs are found on any position of the primers

| SNP position | Loose                                 | Strict                                |
| ------------ | ------------------------------------- | ------------------------------------- |
| -1           | <span style="color:green">PASS</span> | <span style="color:red">FAIL</span>   |
| -2           | <span style="color:red">FAIL</span>   | <span style="color:red">FAIL</span>   |
| -3           | <span style="color:red">FAIL</span>   | <span style="color:red">FAIL</span>   |
| -4           | <span style="color:red">FAIL</span>   | <span style="color:red">FAIL</span>   |
| -5           | <span style="color:green">PASS</span> | <span style="color:red">FAIL</span>   |
| -6 etc       | <span style="color:green">PASS</span> | <span style="color:red">FAIL</span>   |
| 0 found      | <span style="color:green">PASS</span> | <span style="color:green">PASS</span> |

#### Secondary structure filter
##### Template
First the secondary structural elements found in the template on the location of the forward and reverse primer is considered.  If 0 structures are predicted here the primer pair passes for the template secondary structure. If structures are predicted the following citeria are used:
Strict:
- If delta G is lower than -15:
	- If secondary structures are predicted in any of the primers it is very likely that they will occur. => FAIL_template
	- If no secondary structure occur => PASS
- If delta G larger than -15:
	- PASS
Loose:
- If delta G is lower than -15:
	- If secondary structures are predicted in the adapted primer it is very likely that they will occur. If they are in position -1,-2,-3,-4,-5 or -6 => FAIL_template
	- If not on these positions => PASS
- if delta G is higher than -15:
	- PASS


| Altered primer | Common primer | Delta G | Loose                                 | Strict                                |
| -------------- | ------------- | ------- | ------------------------------------- | ------------------------------------- |
| -1 till -6     | 0 predicted   | -10     | <span style="color:green">PASS</span> | <span style="color:green">PASS</span> |
| -1 till -6     | 0 predicted   | -16     | <span style="color:red">FAIL</span>   | <span style="color:red">FAIL</span>   |
| 0 predicted    | -1 till -6    | -10     | <span style="color:green">PASS</span> | <span style="color:green">PASS</span> |
| 0 predicted    | -1 till -6    | -16     | <span style="color:green">PASS</span> | <span style="color:red">FAIL</span>   |
| -1 till -6     | -1 till -6    | -10     | <span style="color:green">PASS</span> | <span style="color:green">PASS</span> |
| -1 till -6     | -1 till -6    | -16     | <span style="color:red">FAIL</span>   | <span style="color:red">FAIL</span>   |
| any            | any           | -10     | <span style="color:green">PASS</span> | <span style="color:green">PASS</span> |
| any            | any           | -16     | depends                               | <span style="color:red">FAIL</span>   |
| 0 predicted    | 0 predicted   | /       | <span style="color:green">PASS</span> | <span style="color:green">PASS</span> |
| -7 till 5'     | any           | <-15    | <span style="color:green">PASS</span> | <span style="color:red">FAIL</span>   |

 ##### Amplicon
 If the amplicon delta G is below -15 it automatically fails regardless (almost never happens; with values this low a structure will always be predicted). When the amplicon detlta G is between -15 and -5 the amplicon only fails when there are secondary structural elements predicted in the amplicon, it passes when there are none.

| Sec str predicted | delta G            | Filter                                |
| ----------------- | ------------------ | ------------------------------------- |
| Yes               | <-15               | <span style="color:red">FAIL</span>   |
| Yes               | -15\< delta G > -5 | <span style="color:red">FAIL</span>   |
| No                | -15\< delta G > -5 | <span style="color:green">PASS</span> |
| Yes               | > -5               | <span style="color:green">PASS</span> |
| No                | > -5               | <span style="color:green">PASS</span> |

#### Validation Filter

Validation will be performed regardless of whether or not this filter is tuned on. While this is useful information to narrow down your options we would recommend looking at this parameter yourself, since this filter is quite strict and will cause most of the lost primer pairs. If the validation gives a specific reason the primer will most likely not be usable. However if no specific reason is given it might be that the primers are usable. Often altered primers get rejected without a clear reason because they are no exact match to the WT template and result in a penalty score that is too high. We would recommend validating this yourself or trying another common primer if very little sequences pass. In summary, we recommend using this filter more to narrow down options than to strictly reject primers. **If a lost of primers are lost because of this step consider turning off this filter.**

The criteria are:
Strict:
- if any of the two primers does not give "considered 1, ok 1":
	- FAIL
Loose:
- if the altered primer does not give "considered 1, ok 1":
	- FAIL
Note: altered primer will often return ok 0 because they penalty score became too high, this does not necessarily mean they are not good if no clear reason is given.


| Altered primer | common primer | Loose                                 | Strict                                | comment                                                             |
| -------------- | ------------- | ------------------------------------- | ------------------------------------- | ------------------------------------------------------------------- |
| ok 0           | ok 1          | <span style="color:red">FAIL</span>   | <span style="color:red">FAIL</span>   | Evaluate primer information                                         |
| ok 1           | ok 0          | <span style="color:green">PASS</span> | <span style="color:red">FAIL</span>   | Evaluate primer information or consider another common primer.      |
| ok 1           | ok 1          | <span style="color:green">PASS</span> | <span style="color:green">PASS</span> | Highest confidence label for this criterium; PASS                   |
| ok 0           | ok 0          | <span style="color:red">FAIL</span>   | <span style="color:red">FAIL</span>   | Consider validation with Primer3plus or using another common primer |
| clear reason   | ok 1          | <span style="color:red">FAIL</span>   | <span style="color:red">FAIL</span>   | FAIL                                                                |
| ok 1           | clear reason  | <span style="color:green">PASS</span> | <span style="color:red">FAIL</span>   | Consider using another common primer                                |
| clear reason   | ok 0          | <span style="color:red">FAIL</span>   | <span style="color:red">FAIL</span>   | FAIL                                                                |
| clear reason   | clear reason  | <span style="color:red">FAIL</span>   | <span style="color:red">FAIL</span>   | Lowest confidence label for this criterium; FAIL                    |



<hr>


## Opening a tsv file in excel

1) Open excel
2) Press Open
<img src="https://github.com/user-attachments/assets/44856f00-9276-4408-a4be-403522a92266" alt="open file" style="width:300px;"/>

4) Browse to the file location (output folder)
<img src="https://github.com/user-attachments/assets/c7ddf87a-38ee-45f4-83be-2eb8430d07ba" alt="browse" style="width:300px"/>

6) Make sure the dropdown in the bottom right is set to: "All Files (\*.\*)"
<img src="https://github.com/user-attachments/assets/a7faa1d9-9406-4a9d-9a60-eea87f743729" alt="all files" style="width:300px"/>

8) Open the file
9) A pop-up will appear, select delimited and next
<img src="https://github.com/user-attachments/assets/28e33377-7908-4522-9eb2-6647462be485" alt="delimiter" style="width:300px"/>

10) Select Tab
<img src="https://github.com/user-attachments/assets/4e9bf009-214c-421d-bae4-35b6f5984338" alt="tab seperated" style="width:300px"/>


11) Next
12) Select general
13) Press Finish
<img src="https://github.com/user-attachments/assets/2548a61a-a8eb-40e8-8b80-f0d9d7546fa2" alt="finish" style="width:300px"/>


#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl=2

/*
====================================================================================================
Pipeline: DMAS
Description: A Nextflow pipeline for the design of double mismatch allele specific primers. 
License: MIT
Copyright (c) 2021 Ghent University
====================================================================================================
*/

/*
====================================================================================================
DEFAULT PARAMETERS (can be overwrittien in config file)
====================================================================================================
*/

// required parameters (input and output)
params.input = "$projectDir/input.txt"
params.outdir = "$projectDir/Output"
// Bowtie index
params.index_bowtie = "$projectDir/Assets/GRCh38/Index_bowtie"
params.index_bowtie_name = "GRCh38_noalt_as"
// Primer3 settings
params.primer_settings = "$projectDir/Assets/GRCh38/Primer3_settings.txt"
// SNP database
params.snp_url = 'http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.bb'
// Primer3 parameters:
params.primer3_diff = 1		// primer min left 3 prime distance
params.primer3_nr = 20		// number of primers to return
// Primer Tm parameters:
params.diff_tm = 2			// melting temp difference between primers
params.min_tm = 58			// min melting temp
params.max_tm = 60			// max melting temp
params.opt_tm = 59			// optimal melting temp
// Primer GC parameters:
params.min_gc = 30			// min GC content
params.opt_gc = 50			// optimal GC content
params.max_gc = 80			// max GC content
// Amplicon length parameters:
params.amp_min = 60			// min amplicon length
params.amp_max = 150 		// max amplicon length
// mispriming library
params.mis_lib = "$projectDir/Assets/humrep_and_simple.txt"		// mispriming library
params.max_mis_lib = 12		// max mispriming library
// filters
params.snp_filter = 'loose'
params.spec_filter = 'loose'
params.sec_str_filter = 'loose'
params.validation_filter = 'loose'
// coordinates included in input file
params.coords = true
params.species = 'human'

// help message
params.help = false

// number of CPUs
params.cpus = 3

// required parameters
input_file = file(params.input)
index_bowtie = file(params.index_bowtie)
primer_settings = file(params.primer_settings)

/*
====================================================================================================
HELP MESSAGE
====================================================================================================
*/

def helpMessage() {
	log.info"""
	Usage:
	
	The typical command for running the pipeline is as follows (standard = default parameters):
	nextflow run DMAS.nf -profile standard --input template.txt --cpus 3
	
	Mandatory nextflow arguments:
	-profile 		set to 'local' when running locally, set to 'singularity' when running on the HPC

	Mandatory pipeline arguments:
	  --input			path to the file (tab separated) with sequence to design primers for (0-based annotation)
	  --index_bowtie		path to bowtie genome index directory
	  --index_bowtie_name	the basename of the Bowtie index to be searched (the name of any of the index files up to but not including the final .1.ebwt / .rev.1.ebwt / ...)
	  --outdir			path to directory where output files should be saved

	Optional pipeline arguments:
	  --coords			true if the input file contains coordinates, false if not (default: true)
	  --cpus			number of CPUs to use in bowtie2 (default: 3) => IMPORTANT: this parameter should be set to the number of CPUs available on the system

	Optional arguments:
	  --primer_settings	path to file with primer3 settings (see primer3 manual)
	  --snp_url			SNP database URL (default: Homo sapiens hg38)
	  --spec_filter		stringency of filtering primer specificity, can be set to 'strict', 'loose' or 'off' (default: loose) => Turn of to run pipeline fast
	  --snp_filter		when set to 'strict', no common SNPs are allowed in primer sequence; when set to 'loose', no SNPs are allowed on position -2,-3,-4, can be turned off with "off"
	  --snp_url			when using a differente species than human, the correct SNP database url should be provided; alternatively, this paramater can be set to 'off' if no SNP database is available
	  --sec_str_filter	when set to 'strict', no secondary structure elements are allowed in primer sequence; when set to 'loose', no secondary structure elements are allowed on position -2,-3,-4,-5,-6 can be turned off with "off"
	  --validation_filter when set to 'strict', both primers need to pass the validation; when set to 'loose', only the specific primer needs to pass the validation; can be turned off with "off"
	  
	  Primer3 settings:
		--primer3_diff	the minimum number of base pairs between the 3' ends of any two left primers (see also primer3 PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE)
		--primer3_nr	number of primers to return (default: 20)
		--min_tm		min melting temp (default: 58)
		--max_tm		max melting temp (default: 60)
		--opt_tm		optimal melting temp (default: 59)
		--diff_tm		melting temp difference between primers (default: 2)
		--min_gc		min GC content (default: 30)
		--max_gc		max GC content (default: 80)
		--opt_gc		optimal GC content (default: 50)
		--amp_min		min amplicon length (default: 50)
		--amp_max		max amplicon length (default: 150)
		--mis_lib		fasta file with mispriming library (default: humrep_and_simple)
		--max_mis_lib	max allowed weighted similarity with any sequence in mispriming library file (default: 12)
	"""
}

// if user defines --help, the error message is printed in the terminal
if (params.help) {
	helpMessage()
	exit 0
}

// Check if parameters are provided correctly
// --input
if (!file(params.input).exists()) {exit 1, "Input file not found: ${params.input}"}
// --outdir
if (!file(params.outdir).exists()) {exit 1, "Output directory not found: ${params.outdir}"}
// --primer_settings
if (!file(params.primer_settings).exists()) {exit 1, "Primer3 settings file not found: ${params.primer_settings}"}
// --index_bowtie
if (!file(index_bowtie).exists()) {exit 1, "Index file not found: ${index_bowtie}"}
// --primer3_diff
if (!params.primer3_diff.toString().isNumber()){exit 1, "Invalid primer3_diff PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE: ${params.primer3_diff}. Valid options: any integer > 0."}
if (params.primer3_diff.toInteger() < 0){exit 1, "Invalid primer3_diff PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE: ${params.primer3_diff}. Valid options: any integer > 0."}
// --primer3_nr
if (!params.primer3_nr.toString().isNumber()){exit 1, "Invalid primer3_nr PRIMER_NUM_RETURN: ${params.primer3_nr}. Valid options: any integer > 0. Caution: setting this parameter to a large value will increase running time."}
if (params.primer3_nr.toInteger() < 0){exit 1, "Invalid primer3_nr PRIMER_NUM_RETURN: ${params.primer3_nr}. Valid options: any integer > 0. Caution: setting this parameter to a large value will increase running time."}
// --min_tm
if (!params.min_tm.toString().isNumber()) {exit 1, " min_tm: ${params.min_tm}. Valid options: any integer > 0."}
// --max_tm
if (!params.max_tm.toString().isNumber()) {exit 1, " max_tm: ${params.max_tm}. Valid options: any integer > 0."}
// --opt_tm
if (!params.opt_tm.toString().isNumber()) {exit 1, " opt_tm: ${params.opt_tm}. Valid options: any integer > 0."}
// --diff_tm
if (!params.diff_tm.toString().isNumber()) {exit 1, " diff_tm: ${params.diff_tm}. Valid options: any integer > 0."}
// --min_gc
if (!params.min_gc.toString().isNumber()) {exit 1, " min_gc: ${params.min_gc}. Valid options: any integer > 0."}
// --max_gc
if (!params.max_gc.toString().isNumber()) {exit 1, " max_gc: ${params.max_gc}. Valid options: any integer > 0."}
// --opt_gc
if (!params.opt_gc.toString().isNumber()) {exit 1, " opt_gc: ${params.opt_gc}. Valid options: any integer > 0."}
// --amp_min
if (!params.amp_min.toString().isNumber()) {exit 1, " amp_min: ${params.amp_min}. Valid options: any integer > 0."}
// --amp_max
if (!params.amp_max.toString().isNumber()) {exit 1, " amp_max: ${params.amp_max}. Valid options: any integer > 0."}
// checking logic
if (params.min_tm.toInteger() > params.max_tm.toInteger() ) {exit 1, " min_tm and max_tm: max_tm (${params.min_tm}) should be > min_tm (${params.max_tm})"}
if (params.opt_tm.toInteger() > params.max_tm.toInteger() || params.opt_tm.toInteger() < params.min_tm.toInteger() ) {exit 1, " opt_tm: ${params.opt_tm} should > min_tm (${params.min_tm}) and < max_tm (${params.max_tm})"}
if (params.min_gc.toInteger() > params.max_gc.toInteger() ) {exit 1, " min_gc and max_gc: max_gc (${params.min_gc}) should be > min_gc(${params.max_gc})"}
if (params.opt_gc.toInteger() > params.max_gc.toInteger() || params.opt_gc.toInteger() < params.min_gc.toInteger() ) {exit 1, " opt_gc: ${params.opt_gc} should > min_gc (${params.min_gc}) and < max_gc (${params.max_gc})"}
// --snp_filter
if (params.snp_filter != "strict" && params.snp_filter != 'loose' && params.snp_filter != 'off'){
	exit 1, "Invalid SNP filter: ${params.snp_filter}. Valid options: 'strict','loose','off'."}
// --spec_filter
if (params.spec_filter != "strict" && params.spec_filter != 'loose' && params.spec_filter != 'off'){
	exit 1, "Invalid specificity filter: ${params.spec_filter}. Valid options: 'strict','loose','off'."}
// --sec_str_filter
if (params.sec_str_filter != "strict" && params.sec_str_filter != 'loose' && params.sec_str_filter != 'off'){
	exit 1, "Invalid secondary structure filter: ${params.sec_str_filter}. Valid options: 'strict','loose','off'."}
// --validation_filter
if (params.validation_filter != "strict" && params.validation_filter != 'loose' && params.validation_filter != 'off'){
	exit 1, "Invalid validation filter: ${params.validation_filter}. Valid options: 'strict','loose','off'."}

log.info """\
==============================================
D M A S   P I P E L I N E
==============================================
OncoRNALab - Arne Blom / Annelien Morlion / Marieke Vromman / Niels Thomas
Github - https://github.com/OncoRNALab/DMAS
Docker - https://hub.docker.com/r/oncornalab/dmas
==============================================
your input file: ${params.input}
your output directory: ${params.outdir}
"""

/*
====================================================================================================
PROCESS 1 - splitting input file
====================================================================================================
*/
// channels
input_file_handle = channel.fromPath(params.input)


// process
process splitInput {
	tag "Splitting input file"
	cpus = params.cpus

	input:
	path('input_file_handle')
	path(bowtie_index)

	output:
	path 'DMAS-*'
	path 'start_time.txt'
	"""
	01_generate_template.py -c $params.coords -i $input_file_handle -p ${task.cpus} -b $bowtie_index/$params.index_bowtie_name
	python3 -c 'from datetime import datetime; print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))' > start_time.txt
	"""
}

/*
====================================================================================================
PROCESS 2 - check coordinates
====================================================================================================
This step is only performed when the species is human!
This tool was developed with the human genome in mind, however if the bowtie index is available for 
other species, the pipeline can be run for other species as well.
====================================================================================================
*/

process checkCoordinates {
	tag "Checking coordinates"
	cpus = params.cpus

	input:
	path(ind_dmas_file_handle)

	output:
	path 'warning.txt'

	when: params.coords == true && params.species == 'human'
	when: params.species == 'human'

	script:
	"""
	02_check_coord.py -i $ind_dmas_file_handle
	"""
}

/*
====================================================================================================
PROCESS 3 - Secondary structure of the template
====================================================================================================
*/

process sec_str_temp {
	tag "Checking secondary structure of the template"

	input:
	path(ind_dmas_file_handle)

	output:
	tuple val("${ind_dmas_file_handle.getBaseName()}"), path('sec-str_*_wt.txt')
	tuple val("${ind_dmas_file_handle.getBaseName()}"), path('sec-str_*_mut.txt')

	script:
	"""
	03_get_sec_str_temp_ViennaRNA.py -t $ind_dmas_file_handle
	"""
}

/*
====================================================================================================
PROCESS 4 - Find common single nucleotide polymorphisms
====================================================================================================
*/

process findSNPs {
	tag "Finding common single nucleotide polymorphisms"

	input:
	path(ind_dmas_file_handle)
	
	output:
	tuple val("${ind_dmas_file_handle.getBaseName()}"), path('snps_*.bed')
	script:
	"""
	04_get_SNPs.py -i $ind_dmas_file_handle -u $params.snp_url
	"""
}

/*
====================================================================================================
THE WORKFLOW
====================================================================================================
*/
workflow {
	// process 1
	splitInput(
		input_file_handle,
		params.index_bowtie)
	// process 2
	checkCoordinates(
		splitInput.out[0].flatten()
	)
	// process 3
	sec_str_temp(
		splitInput.out[0].flatten()
	)
	// process 4
	findSNPs(
		splitInput.out[0].flatten()
	)
}

workflow.onComplete {
	println "\n\t\t\t  Pipeline execution summary\n"+
		"=================================================================================\n\n"+
		"\tPipeline completed at:\t$workflow.complete\n" +
		"\tExecution status:\t${ workflow.success ? 'OK' : 'failed' }\n"+
		"\tNextflow version:\t$nextflow.version\n"+
		"\tExecuted command:\t$workflow.commandLine\n"+
		"\tRun name:\t\t$workflow.runName\n"+
		"\tConfig file:\t\t$workflow.configFiles\n"+
		"\tProfile:\t\t$workflow.profile\n"+
		"\tContainer engine:\t$workflow.containerEngine\n"+
		"\tContainer:\t\t$workflow.container\n"+
		"\tStart time:\t\t$workflow.start\n"+
		"\tCompletion:\t\t$workflow.complete\n"+
		"\tDuration:\t\t$workflow.duration\n"+
		"\tProject directory:\t$workflow.projectDir\n"+
		"\tExit status:\t\t$workflow.exitStatus\n"+
		"================================================================================="
}
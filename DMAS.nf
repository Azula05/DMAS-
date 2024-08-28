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
params.min_left_prime = 1	// min left primer distance
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
params.amp_min = 50			// min amplicon length
params.amp_max = 150 		// max amplicon length
// mispriming library
params.mis_lib = "$projectDir/Assets/GRCh38/humrep_and_simple.txt"		// mispriming library
params.max_mis_lib = 12		// max mispriming library
// salt conentrations
params.dnac = 250
params.na = 50
params.k = 0
params.tris =75 
params.mg = 3 
params.dNTPs = 1.2
params.position = "all"
params.single_MM_Tm = 55
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
	  --primer3_nr		number of primers to return (default: 20)
	  --min_left_prime = 1 min left primer distance
	  --min_tm			min melting temp (default: 58)
	  --max_tm			max melting temp (default: 60)
	  --opt_tm			optimal melting temp (default: 59)
	  --diff_tm			melting temp difference between primers (default: 2)
	  --min_gc			min GC content (default: 30)
	  --max_gc			max GC content (default: 80)
	  --opt_gc			optimal GC content (default: 50)
	  --amp_min			min amplicon length (default: 50)
	  --amp_max			max amplicon length (default: 150)
	  --mis_lib			fasta file with mispriming library (default: humrep_and_simple)
	  --max_mis_lib		max allowed weighted similarity with any sequence in mispriming library file (default: 12)
	
	Temperature prediction:
	  --dnac			DNA concentration nM 	(default: 250)
	  --na				Na+ concentration mM 	(default: 50)
	  --k				K+ concentration mM 	(default: 0)
	  --tris			Tris concentration mM 	(default: 75)
	  --mg				Mg2+ concentration mM	(default: 3)
	  --dNTPs			dNTPs concentration mM	(default: 1.2)
	  --position		2,3,4 or 'all' 			(default: all)
	  --single_MM_Tm	melting temperature for single mismatches °C (default: 55)
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
// --mismatch library
if (!file(params.mis_lib).exists()) {exit 1, "Mispriming library file not found: ${params.mis_lib}"}
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
// --position
if (params.position != "2" && params.position != '3' && params.position != '4' && params.position != 'all'){
	exit 1, "Invalid position: ${params.position}. Valid options: 2,3,4,'all'."
}
// --single_MM_Tm
if (!params.single_MM_Tm.toString().isNumber()) {exit 1, " single_MM_Tm: ${params.single_MM_Tm}. Valid options: any integer > 0."}

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
def input_file_handle = file(params.input)

process splitInput {
	publishDir "$params.outdir/warnings", mode: 'copy', overwrite: false, pattern : 'warning_*.txt'

	tag "Splitting input file"
	cpus = params.cpus

	input:
	path(input_file_handle)
	path(bowtie_index)

	output:
	path('DMAS-*.txt')
	path('start_time.txt')
	path('DMAS_log.txt')

	"""
	01_generate_template.py -c $params.coords -i $input_file_handle -p ${task.cpus} -b $bowtie_index/$params.index_bowtie_name
	python3 -c 'from datetime import datetime; print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))' > start_time.txt
	echo "Name	Total	Specificity_loose	Specificity_strict	SNP_loose	SNP_strict	Sec_str_loose	Sec_str_strict	Sec_str_amp_loose	Sec_str_amp_strict	Validation_loose	Validation_strict	Pass_loose	Pass_strict" > DMAS_log.txt
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
	publishDir "$params.outdir/warnings", mode: 'copy', overwrite: false, pattern : 'warningDMAS-*.txt'

	tag "Checking coordinates"
	cpus = params.cpus

	input:
	path(ind_dmas_file_handle)

	output:
	path 'warningDMAS-*.txt'

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
	03_get_sec_str_temp_ViennaRNA.py -t $ind_dmas_file_handle -T $params.opt_Tm -na $params.na
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
PROCESS 5 - Primer generation
====================================================================================================
*/

process primerGeneration {
	tag "Creating specific primers"

	input:
	path(ind_dmas_file_handle)

	output:
	tuple val("${ind_dmas_file_handle.getBaseName()}"), path('DMAS-*_primers.tsv')

	script:
	"""
	05_primer_generation.py -dnac $params.dnac -Na $params.na -K $params.k -Tris $params.k -Mg $params.mg -dNTPs $params.dNTPs -t $ind_dmas_file_handle -p $params.position -S $params.single_MM_Tm
	"""
}

/*
====================================================================================================
PROCESS 6 - Common primer generation
====================================================================================================
*/
process createTuple {
	tag "Creating tuple"

	input:
	path(ind_dmas_file_handle)

	output:
	tuple val("${ind_dmas_file_handle.getBaseName()}"), path('DMAS-*.tsv')

	script:
	"""
	cat $ind_dmas_file_handle > ${ind_dmas_file_handle.getBaseName()}.tsv
	"""
}

process commonPrimer {
	tag "Creating common primers"

	input:
	path(mis_lib_file)
	tuple val(dmas_id),path(ind_dmas_file_handle), path(ind_dmas_table_handle), path(ind_dmas_snp_handle)
	path(primer3_settings)

	output:
	tuple val(dmas_id), path('Common_FWD.txt')
	tuple val(dmas_id), path('Common_REV.txt')

	script:
	"""
	06_Primer3_common_primer.py -a $params.min_left_prime -b $params.primer3_nr -c $params.min_tm -d $params.max_tm -e $params.opt_tm -f $params.diff_tm -g $params.min_gc -i $params.max_gc -j $params.opt_gc -k $params.amp_min -l $params.amp_max -m $mis_lib_file -M $params.max_mis_lib -dnac $params.dnac -Na $params.na -K $params.k -Tris $params.tris -Mg $params.mg -dNTPs $params.dNTPs -t $ind_dmas_file_handle -p $ind_dmas_table_handle -s $ind_dmas_snp_handle
	primer3_core -p3_settings_file=$primer3_settings < ./Primer3_DMAS-*_common_primer_REV.txt > Common_REV.txt
	primer3_core -p3_settings_file=$primer3_settings < ./Primer3_DMAS-*_common_primer_FWD.txt > Common_FWD.txt
	"""
}

/*
====================================================================================================
PROCESS 7 - Validate primers
====================================================================================================
*/

process primerValidation {
	publishDir "$params.outdir/common_alternatives", mode: 'copy', overwrite: false, pattern : 'DMAS-*_common_alternatives.txt'
	tag "Validating primers"

	input:
	path(primer3_settings)
	tuple val(dmas_id), path(ind_dmas_file_handle), path(ind_dmas_table_handle), path(ind_snp_file), path(temp_sec_str_wt), path(temp_sec_str_mut), path(common_FWD),path(common_REV)

	output:
	tuple val(dmas_id), path('DMAS-*_primer_validated.tsv')
	tuple val(dmas_id), path('DMAS-*_common_alternatives.txt')

	script:
	"""
	07_primer_validation.py -F $common_FWD -R $common_REV -T $ind_dmas_file_handle -s $ind_snp_file -u $temp_sec_str_wt -U $temp_sec_str_mut -o $ind_dmas_table_handle -p $primer3_settings -q yes
	"""
}

/*
====================================================================================================
PROCESS 8 - Get secondary structure of the amplicon
====================================================================================================
*/

process sec_str_amp{

	tag "Checking secondary structure of the amplicon"

	input:
	tuple val(dmas_id), path(ind_dmas_table_handle)

	output:
	tuple val(dmas_id), path('DMAS-*_primers_val_amp.tsv')

	script:
	"""
	08_get_sec_str_amp_ViennaRNA.py -i $ind_dmas_table_handle -T $params.opt_Tm -na $params.na
	"""
}

/*
====================================================================================================
PROCESS 9 - Check the specificity of the primers
====================================================================================================
*/

process spec_primer {
	cpus params.cpus
	maxForks 1

	tag "Checking the specificity of the primers"

	input:
	path(bowtie_index)
	tuple val(dmas_id), path(ind_dmas_table_handle)

	output:
	tuple val(dmas_id), path('DMAS-*_primers_spec.tsv')

	script:
	"""
	09_primer_specificity.py -b $bowtie_index/$params.index_bowtie_name -t $params.cpus -s $params.spec_filter -i $ind_dmas_table_handle
	"""
}

/*
====================================================================================================
PROCESS 10 - Filter the table based on the different filter parameters
====================================================================================================
*/

process filters {
	publishDir "$params.outdir", mode: 'copy', overwrite: true, pattern : 'DMAS_log.txt'
	publishDir "$params.outdir", mode: 'copy', overwrite: true, pattern : 'DMAS-*_primers.tsv'
	publishDir "$params.outdir/Filtered", mode: 'copy', overwrite: true, pattern : 'DMAS-*_filtered_loose.tsv'
	publishDir "$params.outdir/Filtered", mode: 'copy', overwrite: true, pattern : 'DMAS-*_filtered_strict.tsv'

	tag "Filtering primers"

	input:
	tuple val(dmas_id), path(ind_dmas_file_handle), path(ind_dmas_table_handle)
	path(log_file)

	output:
	path('DMAS_log.txt')
	tuple val(dmas_id), path('DMAS-*_filtered_loose.tsv')
	tuple val(dmas_id), path('DMAS-*_filtered_strict.tsv')
	tuple val(dmas_id), path('DMAS-*_primers.tsv')

	script:
	"""
	10_filter.py -i $ind_dmas_file_handle -p $ind_dmas_table_handle -s $params.spec_filter -S $params.snp_filter -t $params.sec_str_filter -v $params.validation_filter -l $log_file
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
	// process 5
	primerGeneration(
		splitInput.out[0].flatten()
	)
	// process 6
	createTuple(
		splitInput.out[0].flatten()
	)
	commonPrimer(
		params.mis_lib,
		createTuple.out[0].join(primerGeneration.out[0]).join(findSNPs.out[0]).groupTuple(),
		params.primer_settings
	)
	// process 7
	primerValidation(
		params.primer_settings,
		createTuple.out[0].join(primerGeneration.out[0]).join(findSNPs.out[0]).join(sec_str_temp.out[0]).join(sec_str_temp.out[1]).join(commonPrimer.out[0]).join(commonPrimer.out[1]).groupTuple()
	)
	// process 8
	sec_str_amp(
		primerValidation.out[0]
	)
	// process 9
	spec_primer(
		params.index_bowtie,
		sec_str_amp.out[0]
	)
	// process 10
	filters(
		createTuple.out[0].join(spec_primer.out[0]).groupTuple(),
		splitInput.out[2]
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
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// required parameters
params.input = "$baseDir/input.txt"
params.outdir = "$baseDir/Output"

// set default parameters
params.settings = "$baseDir/bin/settings_file.txt"
params.index_bowtie = "$baseDir/Data/GRCh38_noalt_as/GRCh38_noalt_as"
params.snp_url = "http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153Common.bb"

params.primer3_diff = 1		// primer min left 3 prime distance
params.primer3_nr = 20		// number of primers to return
params.min_tm = 58			// min melting temp
params.max_tm = 60			// max melting temp
params.opt_tm = 59			// optimal melting temp
params.diff_tm = 2			// melting temp difference between primers
params.min_gc = 30			// min GC content
params.max_gc = 80			// max GC content
params.opt_gc = 50			// optimal GC content
params.amp_min = 60			// min amplicon length
params.amp_max = 150		// max amplicon length
params.mis_lib = "no"		// mispriming library
params.max_mis_lib = 12		// max mispriming library

params.spec_filter = "strict"
params.upfront_filter = "yes"
params.coords = false
params.help = false


// Check if parameters are provided correctly
if (!file(params.input).exists()) {exit 1, "Input file not found: ${params.input}"}
if (!file(params.outdir).exists()) {exit 1, "Output directory not found: ${params.outdir}"}
if (!file(params.settings).exists()) {exit 1, "Primer3 settings file not found: ${params.settings}"}
if (!params.primer3_diff.toString().isNumber()){exit 1, "Invalid primer3_diff PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE: ${params.primer3_diff}. Valid options: any integer > 0."}
if (params.primer3_diff.toInteger() < 0){exit 1, "Invalid primer3_diff PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE: ${params.primer3_diff}. Valid options: any integer > 0."}
if (!params.primer3_nr.toString().isNumber()){exit 1, "Invalid primer3_nr PRIMER_NUM_RETURN: ${params.primer3_nr}. Valid options: any integer > 0. Caution: setting this parameter to a large value will increase running time."}
if (params.primer3_nr.toInteger() < 0){exit 1, "Invalid primer3_nr PRIMER_NUM_RETURN: ${params.primer3_nr}. Valid options: any integer > 0. Caution: setting this parameter to a large value will increase running time."}
if (!params.min_tm.toString().isNumber()) {exit 1, " min_tm: ${params.min_tm}. Valid options: any integer > 0."}
if (!params.max_tm.toString().isNumber()) {exit 1, " max_tm: ${params.max_tm}. Valid options: any integer > 0."}
if (!params.opt_tm.toString().isNumber()) {exit 1, " opt_tm: ${params.opt_tm}. Valid options: any integer > 0."}
if (!params.diff_tm.toString().isNumber()) {exit 1, " diff_tm: ${params.diff_tm}. Valid options: any integer > 0."}
if (!params.min_gc.toString().isNumber()) {exit 1, " min_gc: ${params.min_gc}. Valid options: any integer > 0."}
if (!params.max_gc.toString().isNumber()) {exit 1, " max_gc: ${params.max_gc}. Valid options: any integer > 0."}
if (!params.opt_gc.toString().isNumber()) {exit 1, " opt_gc: ${params.opt_gc}. Valid options: any integer > 0."}
if (!params.amp_min.toString().isNumber()) {exit 1, " amp_min: ${params.amp_min}. Valid options: any integer > 0."}
if (!params.amp_max.toString().isNumber()) {exit 1, " amp_max: ${params.amp_max}. Valid options: any integer > 0."}
if (params.min_tm.toInteger() > params.max_tm.toInteger() ) {exit 1, " min_tm and max_tm: max_tm (${params.min_tm}) should be > min_tm (${params.max_tm})"}
if (params.opt_tm.toInteger() > params.max_tm.toInteger() || params.opt_tm.toInteger() < params.min_tm.toInteger() ) {exit 1, " opt_tm: ${params.opt_tm} should > min_tm (${params.min_tm}) and < max_tm (${params.max_tm})"}
if (params.min_gc.toInteger() > params.max_gc.toInteger() ) {exit 1, " min_gc and max_gc: max_gc (${params.min_gc}) should be > min_gc(${params.max_gc})"}
if (params.opt_gc.toInteger() > params.max_gc.toInteger() || params.opt_gc.toInteger() < params.min_gc.toInteger() ) {exit 1, " opt_gc: ${params.opt_gc} should > min_gc (${params.min_gc}) and < max_gc (${params.max_gc})"}
if (params.spec_filter != "strict" && params.spec_filter != 'loose'){exit 1, "Invalid specificity filter: ${params.spec_filter}. Valid options: 'strict','loose'."}
if (params.upfront_filter != "yes" && params.upfront_filter != "str" && params.upfront_filter != "snp" && params.upfront_filter != "no") {exit 1, "Invalid upfront filter: ${params.upfront_filter}. Valid options: 'yes', 'str', 'snp', 'no'."}


// Help message
def helpMessage() {
	log.info"""
	Example of a typical command for running the pipeline:
	\$ nextflow run dmas.nf --input template.txt

	Required arguments:
	  --input			file with sequence to design primers for
	  --outdir			path to directory where output files should be saved

	Optional options:
	  --coords			provide only if input file contains sequence coordinates to skip the mapping step

	Optional arguments:
	  --settings		path to file with primer3 settings
	  --index_bowtie	path to bowtie2 index files, including the common file names (excluding .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2)
	  --snp_url			SNP database URL (default: Homo sapiens hg38)
	  --spec_filter		stringency of filtering primer specificity, can be set to 'strict' or 'loose' (default: strict)
	  --upfront_filter	When set to 'yes', SNPs and secundary structures are avoided before primer design; when set to 'str', secundary structures are avoided before primer design; when set to 'snp', snp are avoided before primer design; when set to 'no', no filtering before primer design is performed (default: yes)

	  Primer3 settings:
		--primer3_diff	primer min left 3 prime distance (default: 1)
		--primer3_nr	number of primers to return (default: 20)
		--min_tm		min melting temp (default: 58)
		--max_tm		max melting temp (default: 60)
		--opt_tm		optimal melting temp (default: 59)
		--diff_tm		melting temp difference between primers (default: 2)
		--min_gc		min GC content (default: 30)
		--max_gc		max GC content (default: 80)
		--opt_gc		optimal GC content (default: 50)
		--amp_min		min amplicon length (default: 60)
		--amp_max		max amplicon length (default: 150)
		--mis_lib		fasta file with mispriming library (default: empty)
		--max_mis_lib	max allowed weighted similarity with any sequence in mispriming library file (default: 12)
	"""
}

// if user defines --help, the error message is printed in the terminal
if (params.help) {
	helpMessage()
	exit 0
}


log.info """\
==============================================
D M A S   P I P E L I N E
==============================================
output directory : ${params.outdir}
"""


process preprocess_template {
	/* 
	process checks length of input template and throws warning when too short
	process uses Python script to generate fasta file for bowtie2
	*/

	cpus 1
	memory 5.GB

	input:
	path(usr_input)

	output:
	path("seq-*.txt")

	script:
	// -4 to remove [X/Y]
	"""
	num=4
	full_seq_length=\$(head -n 1 ${usr_input} | awk 'NR==1 { print length(\$0);}')
	seq_length="\$((full_seq_length-num))" 
	if (( \$seq_length < 150 )); then
		orange='\033[33m'
		echo "\${orange}WARNING - Your input sequence is shorter than 150 nucleotides. It is recommended to use a sequence of 150 or more nucleotides\${orange}"
	fi

	01_generate_template.py -c $params.coords -i ${usr_input} -b $params.index_bowtie -p ${task.cpus}
	"""
}

process flatten_files{

	input:
	path(sequence)

	output:
	tuple val("${sequence.baseName}"), path(sequence)

	script:
	"echo sequence"
	
}

process get_sec_str_temp {
	/* 
	process uses Python script to retrieve secondary structures in template
	these positions are avoided during primer design
	*/

	input:
	tuple val(sec_str_ID), path(sequence)

	output:
	tuple val(sec_str_ID), path('sec-str_*.txt')

	script:
	"03_get_sec_str_temp.py -t ${sequence}" 
}


process get_SNPs {
	/*
	process uses Python script to retrieve SNP positions from database
	these positions are avoided during primer design
	*/

	input:
	tuple val(SNP_ID), path(sequence)

	output:
	tuple val(SNP_ID), file('snps_seq-*.bed')


	script:
	"05_get_SNPs.py -i ${sequence} -u $params.snp_url"
}


process upfront_filter {
	/*
	process uses Python script to gather SNP and secondary structure 
	information from previous processes and creates an edited primer3 input 
	file to design primers 
	*/

	input:
	tuple val(u_filter_ID), path(sequence), path(snp), path(structure)

	output:
	tuple val(u_filter_ID), path('input-primer3_seq*.txt')


	script:
	"06_upfront_filter.py -a $params.primer3_diff -b $params.primer3_nr -c $params.min_tm -d $params.max_tm -e $params.opt_tm -f $params.diff_tm -g $params.min_gc -i $params.max_gc -j $params.opt_gc -k $params.amp_min -l $params.amp_max -m $params.mis_lib -M $params.max_mis_lib -q $params.upfront_filter -s ${snp} -t ${sequence} -u ${structure}"
}


process primer3_output {
	/*
	process uses primer3 to generate 6 files with primers
	from 6 input files */
	
	publishDir "$params.outdir/primer3_details", mode:"copy",  pattern: "out-input-primer3_seq-*.txt"

	input:
	tuple val(p3_ID), path(in)
	path(settings_ch)

	output:
	path('primer-design-info_seq-*.txt')
	tuple val(p3_ID), path('primer-lists.txt')
	path('input-primer-spec_seq-*.txt')
	tuple val(p3_ID), path('input-sec-str-amplicon*.txt')
	path('out-input-primer3_seq-*.txt')

	script:
	"""
	07_split_and_prepare_primers.py
	cat primer-list_seq-*.txt > primer-lists.txt
	"""
}


process primer_specificity {
	/*
	process will check, using bowtie2, if designed primers are specific
	*/

	// publishDir params.outdir, mode:"copy"
	cpus 1
	memory 5.GB

	input:
	path(in)
	path(in_spec)

	output:
	path('fail-spec_all.txt')


	script:
	"""
	cat ${in_spec} > in_spec_all_primers.txt
	bowtie2 -p ${task.cpus} --xeq --no-hd --no-sq --quiet -x $params.index_bowtie --12 in_spec_all_primers.txt > out-primer-spec_all.sam
	08_primer_spec_filter.py -f $params.spec_filter
	"""
}


process sec_str_amp {
	/*
	process will determine the secondary structures in the amplicons
	*/
	// publishDir params.outdir, mode:"copy"

	input:
	tuple val(sec_str_ID), path(sec_str_amp_input)

	output:
	tuple val(sec_str_ID), path('out-sec-str-amp_seq.txt')


	script:
	"""09_get_sec_str_amp.py
	cat out-sec-str-amp_seq-*.txt > out-sec-str-amp_seq.txt
	"""
}


process filter_primers {
	/*
	process will filter the designed primers on primer specificity and amplicon
	secondary structures
	*/
	publishDir "$params.outdir/primer3_details", mode:"copy", pattern: "filtered-primers_seq-*.txt"

	input:
	tuple val(filter_ID), path(sequence), path(primers), path(sec_str_amp), path(snps), path(sec_str_temp)
	path(specificity)

	output:
	path('filtered-primers_seq-*.txt')
	path('log-file_seq-*.txt')

	script:
	"""
	11_filter.py -p ${primers} -s ${specificity} -S ${snps} -t ${sec_str_amp} -i ${sequence}"""
}


process gather_output {
	/*
	process will create a tab delimited file with best dmas primers
	*/
	publishDir params.outdir, mode:"copy"

	input:
	path(filtered_primers)
	path(log_file)

	output:
	path("dmas_primers.txt")
	path("log_file.txt")

	script:
	"""
	12_gather_output.py -f ${filtered_primers}
	echo -e "seq_ID	design	primer_found	total_primer_pairs	passed	failed_spec	failed_sec_str_amp" | cat - ${log_file} > log_file.txt
	"""
}

usr_input_ch = Channel.fromPath("$params.input").collect()
settings_ch = Channel.fromPath("$params.settings").collect()

// running workflow with processes
workflow{
	preprocess_template(
		usr_input_ch)
	flatten_files(
		preprocess_template.out[0].flatten())
	get_sec_str_temp(
		flatten_files.out)
	get_SNPs(
		flatten_files.out)
	all_info_1 = flatten_files.out.join(get_SNPs.out).join(get_sec_str_temp.out).groupTuple()
	upfront_filter(
		all_info_1)
	primer3_output(
		upfront_filter.out, 
		settings_ch)
	primer_specificity(
		preprocess_template.out[0],
		primer3_output.out[2].collect())
	sec_str_amp(
		primer3_output.out[3])
	all_info_2 = flatten_files.out.join(primer3_output.out[1]).join(sec_str_amp.out).join(get_SNPs.out).join(get_sec_str_temp.out).groupTuple()
	filter_primers(
		all_info_2,
		primer_specificity.out)
	gather_output(
		filter_primers.out[0].collectFile(name: "filtered-primers.txt"),
		filter_primers.out[1].collectFile(name: "log_file_no_header.txt"))
}

#!/usr/bin/env nextflow

/*
 Pipeline to calculate IDR values for biological replicates, generated pseudo-replicates or to check self-consistency of libraries.
 Authors: 
	- Yuvia A. PEREZ RICO <yuvia.perez-rico@embl.de>
*/

log.info "            IDR calculation - version 0.0.2            "
log.info "#######################################################"
log.info "Sample description file	= ${params.sampleInfo}"
log.info "Output directory		= ${params.outDir}"
log.info "File type with peaks		= ${params.peakType}"
log.info "Value use to rank peaks	= ${params.peakRank}"
log.info "Compare biological replicates	= ${params.replicates}"
log.info "Compare pseudo-replicates	= ${params.pseudoReps}"
log.info "Calculate self-consistency	= ${params.selfConsist}"
log.info "\n"

/*
 Validate input parameters
*/

if( !(params.peakType in ['narrowPeak', 'broadPeak', 'bed', 'gff']) ){
	exit 1, "Invalid file type of peaks = ${params.peakType}"
}

if( !(params.peakRank in ['signal.value', 'p.value', 'q.value', 'score']) ){
	exit 1, "Invalid data to rank peaks = ${params.peakRank}"
}

if( !(params.replicates in [0, 1]) ){
	exit 1, "Invalid option = ${params.replicates}"
}

if( !(params.replicates in [0, 1]) ){
	exit 1, "Invalid option = ${params.pseudoReps}"
}

if( !(params.replicates in [0, 1]) ){
	exit 1, "Invalid option = ${params.selfConsist}"
}

/*
 Validate input files
*/

sdFile = file(params.sampleInfo)
if( !sdFile.exists() ){
	exit 1, "The specified sample description file does not exist = ${params.sampleInfo}"
}
log.info "Checking sample description file = $sdFile"

resDir = file(params.outDir)
if( !resDir.exists() && !resDir.mkdirs() ){
	exit 1, "The specified directory to save results cannot be created = ${params.outDir}\n Check file system access permission"
}
log.info "Checking results directory = $resDir"

/*
 Program versions
*/

process get_program_versions{
	publishDir "${params.outDir}/software", mode: 'move'

	output:
	file('programs_version.txt') into programs_version

	"""
	echo nextflow ${nextflow.version} > tmp_version.txt
	samtools --version | grep samtools >> tmp_version.txt
	bedtools --version >> tmp_version.txt
	macs2 --version >> tmp_version.txt
	idr --version | head -1 >> tmp_version.txt
	sort tmp_version.txt > programs_version.txt
	"""

}

/*
 Create channels with the indicated bam or narrowPeak files
*/

if(params.replicates) {
	Channel
		.fromPath(sdFile)
		.splitCsv(header:true)
		.map{ row -> tuple(row.sample1 + "_" + row.sample2, [file(row.nP1), file(row.nP2)]) }
		.set { samples_peaks }
	Channel
		.empty()
		.set { samples_bam }
} else if(params.pseudoReps) {
	Channel
		.fromPath(sdFile)
		.splitCsv(header:true)
		.map{ row -> tuple(row.sample1 + "_" + row.sample2, [file(row.bam1), file(row.bam2)]) }
		.set { samples_bam }
} else if(params.selfConsist) {
	Channel
		.fromPath(sdFile)
		.splitCsv(header:true)
		.map{ row -> tuple(row.sample, file(row.bam)) }
		.set { samples_bam }
}

/*
 Step 1. Merge and split bam files
*/

process split_bams{
	input:
	set val(name), file(bam) from samples_bam
	
	when:
	params.pseudoReps || params.selfConsist

	output:
	set val(name), file('*00.bam') into samples_replicate1
	set val(name), file('*01.bam') into samples_replicate2

	script:
	if(params.pseudoReps) {
		"""
		# Merge
		samtools merge -@ 10 -f ${name}.bam ${bam}
		samtools sort -@ 10 -n -o ${name}_sorted.bam ${name}.bam
		# Split
		samtools view ${name}_sorted.bam | awk 'NR%2 == 1 {f=\$0;next}{print f"separator"\$0}' > ${name}_sorted_joined.sam
		half=\$(echo \$(wc -l ${name}_sorted_joined.sam | cut -f 1 -d " ")/2 | bc)
		shuf ${name}_sorted_joined.sam > ${name}_shuffled_joined.sam
		split --numeric-suffixes=0 -l \$half --filter='sed "s/separator/\\n/" > \$FILE.sam' ${name}_shuffled_joined.sam "${name}_pseudoRep_"
		# Create bam files
		samtools view -H ${name}_sorted.bam > header.txt
		for i in `ls *pseudoRep*.sam`; do cat header.txt \$i | samtools view -@ 10 -b > \${i/sam/bam}; done
		"""
	} else if(params.selfConsist) {
		"""
		# Sort
		samtools sort -@ 10 -n -o ${name}_sorted.bam ${bam}
		# Split
		samtools view ${name}_sorted.bam | awk 'NR%2 == 1 {f=\$0;next}{print f"separator"\$0}' > ${name}_sorted_joined.sam
		half=\$(echo \$(wc -l ${name}_sorted_joined.sam | cut -f 1 -d " ")/2 | bc)
		shuf ${name}_sorted_joined.sam > ${name}_shuffled_joined.sam
		split --numeric-suffixes=0 -l \$half --filter='sed "s/separator/\\n/" > \$FILE.sam' ${name}_shuffled_joined.sam "${name}_selfRep_"
		# Create bam files
		samtools view -H ${bam} > header.txt
		for i in `ls *selfRep*.sam`; do cat header.txt \$i | samtools view -@ 10 -b > \${i/sam/bam}; done
		"""
	}

}

/*
 Step 2. Call peaks - macs2 using fragment information
*/

process macs2_peaks{
	input:
	set val(name), file(bam) from samples_replicate1.mix(samples_replicate2)

	when:
	params.pseudoReps || params.selfConsist

	output:
	set val(name), file("${bam.baseName}_fragments/*peaks.narrowPeak") into samples_partial_peaks

	"""
	macs2 callpeak -f BAMPE -g mm --keep-dup all --outdir ${bam.baseName}_fragments --min-length 100 -n ${bam.baseName} -t ${bam}
	"""

}

if(params.pseudoReps || params.selfConsist) {
	samples_partial_peaks
		.groupTuple()
		.set { samples_peaks }
}

/*
 Step 3. Calculate IDR
*/

process idr{
	publishDir "${params.outDir}/idr/${name}", mode: 'copy'

	input:
	set val(name), file(nPs) from samples_peaks
	
	output:
	file("${name}*") into idr_results

	script:
	if(params.replicates) {
		"""
		idr --samples ${nPs} --input-file-type ${params.peakType} --rank ${params.peakRank} --output-file ${name}_idr --plot --log-output-file ${name}.idr.log
		"""
	} else if(params.pseudoReps || params.selfConsist) {
		"""
		idr --samples ${nPs} --input-file-type ${params.peakType} --rank ${params.peakRank} --output-file ${name}_pseudoReps_idr \
			--plot --log-output-file ${name}_pseudoReps.idr.log
		"""
	}

}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Main results are saved in ${params.outDir}\n" : "There was an error during the execution, check log files." )
}



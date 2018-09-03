#!/bin/bash

ulimit -n 2048
SHOW_HELP=0
while getopts ":i:g:r:o:b:u:p:n:h:" opt; do
  case $opt in
    i) 	RRFILE="$OPTARG"
    ;;
    g) 	GFFFILE="$OPTARG"
    ;;
	r) 	REFDIR="$OPTARG"
	;;
    o) 	OUTPATH="$OPTARG"
    ;;
    b) 	BCFILE="$OPTARG"
    ;;
    u) 	UBDPATH="$OPTARG"
    ;;
    p) 	BUCKETPATH="$OPTARG"
    ;;
    n) 	NCORES="$OPTARG"
    ;;
	h) 	SHOW_HELP=1
	;;
    \?) SHOW_HELP=1
    ;;
  esac
done

help_readme() {
#Help
	echo "

##############################################################

simpleST version 1.0
Copyright (c) 2018 Martin Fahrenberger.
This wrapper is free software and comes with ABSOLUTELY NO WARRANTY.

##############################################################

This pipeline processes .sra data files as produced using the Spatial Transciptomics protocol by Stahl and Salmen et al. (STAHL, Patrik L., et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. Science, 2016).

##############################################################

This pipeline uses the following software:

fastq-dump:
Part of NCBI's SRA-toolkit

UBD:
Barcode demultiplexer (COSTEA, Paul Igor; LUNDEBERG, Joakim; AKAN, Pelin. TagGD: fast and accurate software for DNA Tag generation and demultiplexing. PLoS One, 2013)

barcode_buckets_v2.py:
custom Python script for efficient sorting of reads by barcode

STAR:
The STAR aligner (DOBIN, Alexander, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 2013)

featureCounts:
Read counter (LIAO, Yang; SMYTH, Gordon K.; SHI, Wei. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 2013)

samtools:
general SAM and BAM file handling (LI, Heng, et al. The sequence alignment/map format and SAMtools. Bioinformatics, 2009)

umi-tools:
UMI deduplication (SMITH, Tom Sean; HEGER, Andreas; SUDBERY, Ian. UMI-tools: Modelling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome research, 2017)

##############################################################
USAGE EXAMPLE:
##############################################################

simpleST.sh -i <SRA-Input-File> -g <reference_gff> -r <reference_STAR_genome_index_dir> -o <outout_dir> -b <barcode_file> -u <UBD_directory> -p <barcode_buckets_path> -n <number_of_cores>

-i		SRA Input file from Spatial_Transcriptomics experiment

-g		Reference genome gff file (same as used for genome index build)

-r		STAR genome index as built with STAR --runMode genomeGenerate

-o		Output Directory

-b		Barcode File containing the coordinates for each barcodes

-u		Path to UBD Installation

-p		Path to python script barcode_buckets_v2.py

-n		Number of cores to use

-h		Shows this help message
"

}
Run_pipeline(){

#Execute the pipeline

	mkdir -p $OUTPATH/fastq

	RRNAME=$(echo $RRFILE | rev | cut -f1 -d/ | cut -f2 -d. | rev)

	echo processing $RRNAME
	mkdir -p $OUTPATH/fastq/$RRNAME
	echo "$(date) Starting fastq-dump"
	fastq-dump --unaligned --split-3 $RRFILE -O $OUTPATH/fastq/$RRNAME/

	cd $OUTPATH/fastq/$RRNAME/
		echo "$(date) Starting split"
		split -l 10000000 -d --additional-suffix=.fastq  ${RRNAME}_1.fastq ${RRNAME}_1_
		split -l 10000000 -d --additional-suffix=.fastq  ${RRNAME}_2.fastq ${RRNAME}_2_
		rm ${RRNAME}_1.fastq
		rm ${RRNAME}_2.fastq
		FILESTMP=($(ls ./*.fastq | grep _1_))
		echo "$(date) Starting UBD demulitplexing"
		for k in ${FILESTMP[@]}; do
		    $UBDPATH/findIndexes -p ${k//_1_/_2_} $BCFILE ./$k ./${k//_1_/}_bc -l 18 >> ./UBD_mapping_stats.txt 
		    rm $k
		    rm ${k//_1_/_2_}
		done
		echo "$(date) Starting python sorting script"
		python3 $BUCKETPATH ./
		rm *_bc*
	cd ../../

	echo "$(date) Starting fastx_trimmer"

	BARCODE_FILES=($(find $OUTPATH/fastq/$RRNAME/assigned/ -name "*_1_.fastq"))

	for i in ${BARCODE_FILES[@]}; do
		fastx_trimmer -f 19 -l 26 -i $i  -o $i.trimmed
		rm $i
	done
	echo "$(date) Starting umi_tools extract"
	TRIMMED_FILES=($(find $OUTPATH/fastq/$RRNAME/assigned/ -name "*_1_.fastq.trimmed"))
	for i in ${TRIMMED_FILES[@]]}; do
		TMP_FILE_2=${i//_1_/_2_}
		umi_tools extract --bc-pattern=NNNNNNNN \
		              --stdin $i \
		              --stdout $i.extracted \
		              --read2-in ${TMP_FILE_2//.trimmed/} \
		              --read2-out=${TMP_FILE_2//.trimmed/}.extracted \
		              >> $OUTPATH/fastq/$RRNAME/assigned/extraction_stats.txt
		rm $i.extracted
		rm $i
		rm ${TMP_FILE_2//.trimmed/}
	done


	echo "$(date) Starting STAR alignments"

	EXTRACTED_FILES=($(find $OUTPATH/fastq/$RRNAME/assigned/ -name "*_2_.fastq.extracted"))
	EXTRACTED_NAMES=($(find $OUTPATH/fastq/$RRNAME/assigned/ -name "*_2_.fastq.extracted" | rev | cut -f1 -d/ | rev | cut -f1 -d_))


	mkdir -p $OUTPATH/alignment/STAR/$RRNAME/
	for j in `seq 1 ${#EXTRACTED_FILES[@]}`; do
		i=$(expr $j - 1)
		mkdir -p $OUTPATH/alignment/STAR/$RRNAME/${EXTRACTED_NAMES[$i]}/
		STAR --genomeDir $REFDIR --readFilesIn ${EXTRACTED_FILES[$i]} --runThreadN $NCORES --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $OUTPATH/alignment/STAR/$RRNAME/${EXTRACTED_NAMES[$i]}/
	done



	echo "$(date) Starting umi_tools dedup (with featureCounts and samtools)"

	ASSIGNED_FOLDERS=($(ls -d $OUTPATH/alignment/STAR/$RRNAME/*))
	for i in ${ASSIGNED_FOLDERS[@]}; do
		featureCounts -T $NCORES -R BAM -g gene -a $GFFFILE -o $i/gene_assigned.txt $i/Aligned.sortedByCoord.out.bam
		samtools sort $i/Aligned.sortedByCoord.out.bam.featureCounts.bam $i/assigned_sorted
		samtools index $i/assigned_sorted.bam
		umi_tools dedup -I $i/assigned_sorted.bam --output-stats=$i/deduplicated -S $i/deduplicated.bam
		featureCounts -T $NCORES -t exon -g gene -a $GFFFILE -o $i/counts.txt $i/deduplicated.bam
	done
}

if [ "$SHOW_HELP" -eq 1 ];
	then
	help_readme
elif [ $# -lt 8 ];
	then
	help_readme
else
	Run_pipeline
fi

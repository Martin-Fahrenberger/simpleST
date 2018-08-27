# simpleST

simpleST version 1.0
Copyright (c) 2018 Martin Fahrenberger.
This wrapper is free software and comes with ABSOLUTELY NO WARRANTY.

This pipeline processes .sra data files as produced using the Spatial Transciptomics protocol by Stahl and Salmen et al. (STAHL, Patrik L., et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. Science, 2016).

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

USAGE EXAMPLE:

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

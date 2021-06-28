#!/bin/bash
set -euo pipefail

if [ ! -d "qc/fastq-output" ]; then
	echo "running fastq analysis"
	mkdir -p qc/fastq-output/
	fastqc fastq_data/SRR1518011_1.fastq.gz fastq_data/SRR1518011_2.fastq.gz -o qc/fastq-output/
	#unzip qc/fastq-output/SRR1518011_1_fastqc.zip -d qc/fastq-output/zip_output/
	#unzip qc/fastq-output/SRR1518011_2_fastqc.zip -d qc/fastq-output/zip_output/
else
	echo "fastq analysis already run"
fi

#if [ ! -d "qc/trim" ]; then
#	echo "trimming adapters and primers"
#	mkdir -p qc/trim/logs
#	mkdir -p qc/trim/paired
#	mkdir -p qc/trim/unpaired
#	TrimmomaticPE -threads 1 -phred33 -trimlog qc/trim/logs/trimm_logfile \
#	fastq_data/SRR61518011_1.fastq.gz fastq_data/SRR61518011_2.fastq.gz \
#	qc/trim/paired/SRR61518011_1_paired.fastq.gz qc/trim/unpaired/SRR61518011_2_unpaired.fastq.gz \
#	qc/trim/paired/SRR61518011_2_paired.fastq.gz qc/trim/unpaired/SRR61518011_2_unpaired.fastq.gz \
#	ILLUMINACLIP:../resources/primers_adapters.fa:2:30:10 MINLEN:36
#else
#	echo "trimming adapters and primers has already occured"
#fi
#
#if [ ! -d "qc/trim2" ]; then
#	echo "trimming low quality reads"
#	mkdir -p qc/trim2/logs
#	mkdir -p qc/trim2/paired
#	mkdir -p qc/trim2/unpaired
#	TrimmomaticPE -threads 1 -phred33 -trimlog qc/trim2/logs/trimm_logfile \
#	qc/trim/paired/SRR61518011_1_paired.fastq.gz qc/trim/paired/SRR61518011_2_paired.fastq.gz \
#	qc/trim2/paired/SRR61518011_1_PE_trimmed_adapter_removed.fastq.gz qc/trim2/unpaired/SRR61518011_1_UP_trimmed_adapter_removed.fastq.gz \
#	qc/trim2/paired/SRR61518011_2_PE_trimmed_adapter_removed.fastq.gz qc/trim2/unpaired/SRR61518011_2_UP_trimmed_adapter_removed.fastq.gz \
#	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#else
#        echo "trimming low quality reads has already occured"
#fi

#if [ ! -d "qc/kmer_freq" ]; then
#	echo "calculating kmer frequency"
#	mkdir -p qc/kmer_freq
#	mv fastq-files.txt qc/kmer_freq
#	cd qc/kmer_freq	
#	KmerFreq_AR -k 16 -t 1 -q 33 -p Error_Corr fastq-files.txt > kmerfreq16.log 2> kmerfreq16.err
#	cd ../../
#else
#	echo "kmer frequency already calculated"
#fi
#
#if [ ! -d "qc/kmer_freq/corrector" ]; then
#	echo "performing correction"
#	mkdir -p qc/kmer_freq/corrector
#	cd qc/kmer_freq
#	Corrector_AR -k 16 -Q 33 -t 1 -o 1 Error_Corr.freq.cz Error_Corr.freq.cz.len fastq-files.txt > Corr16.log 2> Corr16.err
#	mv Corr16* corrector/
#	cd ../../
#else
#	echo "correction already done"
#fi

#!/bin/bash
set -euo pipefail

if [ ! -f "../resources/reference/reference-genome.fasta.bwt" ]; then
	echo "indexing reference genome"
	bwa index /home/jovyan/resources/reference/reference-genome.fasta
else
	echo "genome already indexed"
fi

if [ ! -d "alignment/bwa_mem" ]; then
	echo "aligning sequence"
	mkdir -p alignment/bwa_mem/
	cd alignment/bwa_mem/
	bwa mem -t 8 -R '@RG\tID:identifier\tLB:library\tPL:platform\tPU:platform-unit\tSM:sample-name' \
	../../../resources/reference/reference-genome.fasta \
	../../fastq_data/SRR1518011_1.fastq.gz \
	../../fastq_data/SRR1518011_2.fastq.gz > initial-output.sam
	cd ../../
else
	echo "alignment already done"
fi

if [ ! -f "alignment/bwa_mem/initial-output.bam" ]; then
	echo "converting sam to bam"
	cd alignment/bwa_mem/
	samtools view -b initial-output.sam -o initial-output.bam
	cd ../../
else
	echo "conversion already happened"
fi

if [ ! -d "alignment/sorted_files/" ]; then
	echo "sorting bam"
	mkdir -p alignment/sorted_files/
	cd alignment/sorted_files/
	samtools sort ../bwa_mem/initial-output.bam -o sorted-file.bam
	cd ../../
else
	echo "bam already sorted"
fi

if [ ! -f "alignment/sorted_files/sorted-file.bam.bai" ]; then
	echo "indexing bam"
	cd alignment/sorted_files/
	samtools index sorted-file.bam
	cd ../../
else
	echo "bam already indexed"
fi

if [ ! -f "alignment/flagstats/flagstat.txt" ]; then
	echo "generating stats"
	mkdir -p alignment/flagstats/
	cd alignment/flagstats/
	samtools flagstat ../sorted_files/sorted-file.bam > flagstat.txt
	cd ../../
else
	echo "stats already generated"
fi

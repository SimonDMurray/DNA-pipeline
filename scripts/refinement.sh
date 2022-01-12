#!/bin/bash
set -euo pipefail

if [ ! -d "refinement/duplicates" ]; then
	echo "marking duplicates"
	mkdir -p refinement/duplicates
	cd refinement/duplicates
	java -jar ~/bin/picard/build/libs/picard.jar MarkDuplicates INPUT=../../alignment/sorted_files/sorted-file.bam OUTPUT=removed-duplicates.bam METRICS_FILE=duplicate-metrics.txt
	cd ../../
else
	echo "duplicates marked"
fi

if [ ! -f "refinement/duplicates/removed-duplicates.bam.bai" ]; then
	echo "indexing removed duplicates bam file"
	cd refinement/duplicates
	samtools index removed-duplicates.bam
	cd ../../
else
	echo "removed duplicates bam already indexed"
fi


if [ ! -d "refinement/quality-map-report" ]; then
	echo "producing qualimap qc report"
	cd refinement
	qualimap bamqc -bam duplicates/removed-duplicates.bam -outdir quality-map-report
	cd ../
else
	echo "qualimap qc report already generated"
fi

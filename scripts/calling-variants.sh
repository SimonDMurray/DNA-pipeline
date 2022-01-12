#!/bin/bash
set -euo pipefail

if [ ! -f "../resources/reference/reference-genome.fasta.fai" ]; then
	echo "indexing reference genome"
	samtools faidx ../resources/reference/reference-genome.fasta
else
	echo "genome already indexed"
fi

if [ ! -f "../resources/reference/reference-genome.dict" ]; then
	echo "creating reference dictionary"
	java -jar bin/picard/build/libs/picard.jar CreateSequenceDictionary -R ../resources/reference/reference-genome.fasta
else
	echo "dictionary already made"
fi

if [ ! -d "calling_variants/gatk-dir/" ]; then
	echo "recalibrating data"
	mkdir -p calling_variants/gatk-dir/
	cd calling_variants/gatk-dir/
	gatk BaseRecalibrator -I ../../refinement/duplicates/removed-duplicates.bam -R ../../../resources/reference/reference-genome.fasta \
	--known-sites ../../../resources/resources-broad-hg38-v0-Homo_sapiens_assembly38.known_indels.vcf \
	--known-sites ../../../resources/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf -O recal_data.table 
	cd ../../
else
	echo "data already recalibrated"
fi

if [ ! -f "calling_variants/gatk-dir/rescored.bam" ]; then
	echo "apply base quality score recalibration"
	cd calling_variants/gatk-dir/
	gatk ApplyBQSR -R ../../../resources/reference/reference-genome.fasta -I ../../refinement/duplicates/removed-duplicates.bam -bqsr recal_data.table -O rescored.bam 
	cd ../../
else
	echo "bqsr already applied"
fi

if [ ! -f "calling_variants/gatk-dir/gatk-file.vcf" ]; then
	echo "haplotype caller running"
	cd calling_variants/gatk-dir/
	gatk HaplotypeCaller -R ../../../resources/reference/reference-genome.fasta -I rescored.bam -mbq 20 --minimum-mapping-quality 50 -O gatk-file.vcf
	cd ../../
else
	echo "halpotype already ran"
fi

if [ ! -f "calling_variants/gatk-dir/selected-gatk.vcf" ]; then
	echo "selecting variants"
	cd calling_variants/gatk-dir/
	gatk SelectVariants -R ../../../resources/reference/reference-genome.fasta --variant gatk-file.vcf -O initial-selected-gatk.vcf --select-type SNP
	grep -w '^#\|^#CHROM\|chr[1-9]\|chr[1-2][0-9]\|chr[X,Y,M]' initial-selected-gatk.vcf > selected-gatk.vcf	
	cd ../../
else
	echo "variants already selected"
fi

if [ ! -f "calling_variants/freebayes-dir/freebayes-file.vcf" ]; then
	echo "running freebayes"
	mkdir -p calling_variants/freebayes-dir/
	cd calling_variants/freebayes-dir/
	freebayes -q 20 -m 50 -u -f ../../../resources/reference/reference-genome.fasta ../gatk-dir/rescored.bam > freebayes-initial-file.vcf
	grep -w '^#\|^#CHROM\|chr[1-9]\|chr[1-2][0-9]\|chr[X,Y,M]' freebayes-initial-file.vcf > freebayes-file.vcf
	cd ../../
else
	echo "freebayes already ran"
fi

if [ ! -f "calling_variants/compared/compare-vcf-files.diff.sites_in_files" ]; then
	echo "comparing variant calling"
	mkdir -p calling_variants/compared/
	cd calling_variants/compared/
	vcftools --vcf ../gatk-dir/selected-gatk.vcf --diff ../freebayes-dir/freebayes-file.vcf --diff-site --out compare-vcf-files 
	cd ../../
else
	echo "variant calling already compared"
fi

if [ ! -f "calling_variants/filtering/filtered-vcf.recode.vcf" ]; then
        echo "filter freebayes vcf file"
        mkdir -p calling_variants/filtering/
        cd calling_variants/filtering/
        vcftools --vcf ../gatk-dir/selected-gatk.vcf --minDP 3 --minQ 20 --out temp --recode --recode-INFO-all
        vcftools --vcf temp.recode.vcf --max-missing 1 --out filtered-vcf --recode --recode-INFO-all
        cd ../../
else
        echo "vcf filtering already done"
fi

#!/bin/bash
set -euo pipefail

if [ ! -d "../resources/snpeff-genome" ]; then
	echo "downloading snpeff genome"
	mkdir -p ../resources/snpeff-genome
	java -jar ~/bin/snpEff/snpEff.jar download GRCh38.99
	cp -r ~/bin/snpEff/data/GRCh38.99 ../resources/snpeff-genome/
else
	echo "genome already downloaded"
fi

if [ ! -d "annotations" ]; then
	echo "annotating vcf file"
	mkdir -p annotations
	cd annotations
	java -jar ~/bin/snpEff/snpEff.jar -c ~/bin/snpEff/snpEff.config GRCh38.99 ../calling_variants/filtering/filtered-vcf.recode.vcf > snpeff.vcf
      	cd ../
else
 	echo "vcf already annotated"
fi

if [ ! -d "annotations/CYP2C19" ]; then
	echo "selecting CYP2C19 info"
	mkdir -p annotations/CYP2C19
	cd annotations/CYP2C19
	grep CYP2C19 ../snpeff.vcf > CYP2C19.vcf
	cd ../../
else
	echo "CYP2C19 already selected"
fi

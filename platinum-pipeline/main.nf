#!/home/jovyan/bin/ nextflow

process getFastq {
	
	echo true

	shell:
	'''
	mkdir -p fastq_data
	if [ ! -f "fastq_data/SRR1518011_1.fastq.gz" ]; then
        	echo "fastq_data/SRR1518011_1.fastq.gz does not exist. Downloading now."
        	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/001/SRR1518011/SRR1518011_1.fastq.gz -O fastq_data/SRR1518011_1.fastq.gz
	else
        	echo "fastq_data/SRR1518011_1.fastq.gz exists. Skipping download."
	fi
	if [ ! -f "fastq_data/SRR1518011_2.fastq.gz" ]; then
        	echo "fastq_data/SRR1518011_2.fastq.gz does not exist. Downloading now."
        	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/001/SRR1518011/SRR1518011_2.fastq.gz -O fastq_data/SRR1518011_2.fastq.gz
	else
        	echo "fastq_data/SRR1518011_2.fastq.gz exists. Skipping download."
	fi
	'''
}

process runFastqc {
	
	echo true
	
	shell:
	'''
	PATH="/usr/lib/:/opt/conda/condabin:/usr/local/go/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/lib/rstudio-server/bin"
	if [ ! -d "qc/fastq-output" ]; then
        	echo "running fastq analysis"
        	mkdir -p qc/fastq-output/
        	fastqc fastq_data/SRR1518011_1.fastq.gz fastq_data/SRR1518011_2.fastq.gz -o qc/fastq-output/
        	unzip qc/fastq-output/SRR1518011_1_fastqc.zip -d qc/fastq-output/zip_output/
        	unzip qc/fastq-output/SRR1518011_2_fastqc.zip -d qc/fastq-output/zip_output/
	else
        	echo "fastq analysis already run"
	fi
	PATH="/opt/conda/bin:/opt/conda/condabin:/usr/local/go/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/lib/rstudio-server/bin"
	'''
}

/*
process removePrimersAdapters {
	
	echo true
	
	shell:
	'''
	if [ ! -d "qc/trim" ]; then
      		echo "trimming adapters and primers"
      	 	mkdir -p qc/trim/logs
       		mkdir -p qc/trim/paired
       		mkdir -p qc/trim/unpaired
       		TrimmomaticPE -threads 1 -phred33 -trimlog qc/trim/logs/trimm_logfile \
       		fastq_data/SRR61518011_1.fastq.gz fastq_data/SRR61518011_2.fastq.gz \
       		qc/trim/paired/SRR61518011_1_paired.fastq.gz qc/trim/unpaired/SRR61518011_2_unpaired.fastq.gz \
       		qc/trim/paired/SRR61518011_2_paired.fastq.gz qc/trim/unpaired/SRR61518011_2_unpaired.fastq.gz \
       		ILLUMINACLIP:../resources/primers_adapters.fa:2:30:10 MINLEN:36
	else
       		echo "trimming adapters and primers has already occured"
	fi
	'''
}

process trimLowQuality {

	echo true

	shell:
	'''
	if [ ! -d "qc/trim2" ]; then
       		echo "trimming low quality reads"
       		mkdir -p qc/trim2/logs
       		mkdir -p qc/trim2/paired
       		mkdir -p qc/trim2/unpaired
       		TrimmomaticPE -threads 1 -phred33 -trimlog qc/trim2/logs/trimm_logfile \
       		qc/trim/paired/SRR61518011_1_paired.fastq.gz qc/trim/paired/SRR61518011_2_paired.fastq.gz \
       		qc/trim2/paired/SRR61518011_1_PE_trimmed_adapter_removed.fastq.gz qc/trim2/unpaired/SRR61518011_1_UP_trimmed_adapter_removed.fastq.gz \
       		qc/trim2/paired/SRR61518011_2_PE_trimmed_adapter_removed.fastq.gz qc/trim2/unpaired/SRR61518011_2_UP_trimmed_adapter_removed.fastq.gz \
       		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	else
        	echo "trimming low quality reads has already occured"
	fi
	'''
}
*/

process indexReference1 {
	
	echo

	shell:
	'''
	if [ ! -f "../resources/reference/reference-genome.fasta.bwt" ]; then
        	echo "indexing reference genome"
        	bwa index /home/jovyan/resources/reference/reference-genome.fasta
	else
        	echo "genome already indexed"
	fi
	'''
}

process alignSequence {

	echo true

	shell:
	'''
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
	'''
}

process convertSam {

	echo true

	shell:
	'''
	if [ ! -f "alignment/bwa_mem/initial-output.bam" ]; then
        	echo "converting sam to bam"
        	cd alignment/bwa_mem/
        	samtools view -b initial-output.sam -o initial-output.bam
        	cd ../../
	else
        	echo "conversion already happened"
	fi
	'''
}

process alignBam {

	echo true

	shell:
	if [ ! -d "alignment/sorted_files/" ]; then
        	echo "sorting bam"
        	mkdir -p alignment/sorted_files/
        	cd alignment/sorted_files/
        	samtools sort ../bwa_mem/initial-output.bam -o sorted-file.bam
        	cd ../../
	else
        	echo "bam already sorted"
	fi
	'''
}

process indexBam1 {

	echo true

	shell:
	'''
	if [ ! -f "alignment/sorted_files/sorted-file.bam.bai" ]; then
        	echo "indexing bam"
        	cd alignment/sorted_files/
        	samtools index sorted-file.bam
        	cd ../../
	else
        	echo "bam already indexed"
	fi
	'''
}

process runFlagstat {

	echo true

	shell:
	'''
	if [ ! -f "alignment/flagstats/flagstat.txt" ]; then
        	echo "generating stats"
        	mkdir -p alignment/flagstats/
        	cd alignment/flagstats/
        	samtools flagstat ../sorted_files/sorted-file.bam > flagstat.txt
        	cd ../../
	else
        	echo "stats already generated"
	fi
	'''
}

process markingDuplicates {

	echo true

	shell:
	'''
	if [ ! -d "refinement/duplicates" ]; then
        	echo "marking duplicates"
        	mkdir -p refinement/duplicates
        	cd refinement/duplicates
        	java -jar ~/bin/picard/build/libs/picard.jar MarkDuplicates INPUT=../../alignment/sorted_files/sorted-file.bam OUTPUT=removed-duplicates.bam METRICS_FILE=duplicate-metrics.txt
        	cd ../../
	else
        	echo "duplicates marked"
	fi
	'''
}

process indexBam2 {
	
	echo true

	shell:
	'''
	if [ ! -f "refinement/duplicates/removed-duplicates.bam.bai" ]; then
        	echo "indexing removed duplicates bam file"
        	cd refinement/duplicates
        	samtools index removed-duplicates.bam
        	cd ../../
	else
        	echo "removed duplicates bam already indexed"
	fi
	'''
}

process runQualimap {

	echo true

	shell:
	'''
	if [ ! -d "refinement/quality-map-report" ]; then
        	echo "producing qualimap qc report"
        	cd refinement
        	qualimap bamqc -bam duplicates/removed-duplicates.bam -outdir quality-map-report
        	cd ../
	else
        	echo "qualimap qc report already generated"
	fi
	'''
}


process indexReference2 {

	echo true

	shell:
	'''
	if [ ! -f "../resources/reference/reference-genome.fasta.fai" ]; then
        	echo "indexing reference genome"
        	samtools faidx ../resources/reference/reference-genome.fasta
	else
        	echo "genome already indexed"
	fi
	'''
}

process createReferenceDict {

	echo true

	shell:
	'''
	if [ ! -f "../resources/reference/reference-genome.dict" ]; then
        	echo "creating reference dictionary"
        	java -jar bin/picard/build/libs/picard.jar CreateSequenceDictionary -R ../resources/reference/reference-genome.fasta
	else
        	echo "dictionary already made"
	fi
	'''
}

process recalibrateData {

	echo true

	shell:
	'''
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
	'''
}

process rescoreData {

	echo true

	shell:
	'''
	if [ ! -f "calling_variants/gatk-dir/rescored.bam" ]; then
        	echo "apply base quality score recalibration"
        	cd calling_variants/gatk-dir/
        	gatk ApplyBQSR -R ../../../resources/reference/reference-genome.fasta -I ../../refinement/duplicates/removed-duplicates.bam -bqsr recal_data.table -O rescored.bam
        	cd ../../
	else
        	echo "bqsr already applied"
	fi
	'''
}

process haplotypeCalling {

	echo true

	shell:
	'''
	if [ ! -f "calling_variants/gatk-dir/gatk-file.vcf" ]; then
        	echo "haplotype caller running"
        	cd calling_variants/gatk-dir/
        	gatk HaplotypeCaller -R ../../../resources/reference/reference-genome.fasta -I rescored.bam -mbq 20 --minimum-mapping-quality 50 -O gatk-file.vcf
        	cd ../../
	else
        	echo "halpotype already ran"
	fi
	'''
}

process selectingVariants {

	echo true

	shell:
	if [ ! -f "calling_variants/gatk-dir/selected-gatk.vcf" ]; then
        	echo "selecting variants"
        	cd calling_variants/gatk-dir/
        	gatk SelectVariants -R ../../../resources/reference/reference-genome.fasta --variant gatk-file.vcf -O initial-selected-gatk.vcf --select-type SNP
        	grep -w '^#\|^#CHROM\|chr[1-9]\|chr[1-2][0-9]\|chr[X,Y,M]' initial-selected-gatk.vcf > selected-gatk.vcf
        	cd ../../
	else
        	echo "variants already selected"
	fi
	'''
}

process runFreebayes {

	echo true

	shell:
	'''
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
	'''
}

process compareCalling {

	echo true

	shell:
	if [ ! -f "calling_variants/compared/compare-vcf-files.diff.sites_in_files" ]; then
        	echo "comparing variant calling"
        	mkdir -p calling_variants/compared/
        	cd calling_variants/compared/
        	vcftools --vcf ../gatk-dir/selected-gatk.vcf --diff ../freebayes-dir/freebayes-file.vcf --diff-site --out compare-vcf-files \
        	cd ../../
	else
        	echo "variant calling already compared"
	fi
	'''
}

process filterVariants {

	echo true

	shell:
	'''
	if [ ! -f "calling_variants/filtering/filtered-vcf.recode.vcf" ]; then
        	echo "filter freebayes vcf file"
        	mkdir -p calling_variants/filtering/
        	cd calling_variants/filtering/
        	vcftools --vcf ../freebayes-dir/freebayes-file.vcf --minDP 3 --minQ 20 --out temp --recode --recode-INFO-all
        	vcftools --vcf temp.recode.vcf --max-missing 1 --out filtered-vcf --recode --recode-INFO-all
        	cd ../../
	else
        	echo "vcf filtering already done"
	fi
}

process downloadSnpeffGenome {

	echo true

	shell:
	if [ ! -d "../resources/snpeff-genome" ]; then
	        echo "downloading snpeff genome"
        	mkdir -p ../resources/snpeff-genome
        	java -jar ~/bin/snpEff/snpEff.jar download GRCh38.99
       		cp -r ~/bin/snpEff/data/GRCh38.99 ../resources/snpeff-genome/
	else
        	echo "genome already downloaded"
	fi
	'''
}

process annotateVCF {

	echo true

	shell:
	if [ ! -d "annotations" ]; then
        	echo "annotating vcf file"
        	mkdir -p annotations
        	cd annotations
        	java -jar ~/bin/snpEff/snpEff.jar -c ~/bin/snpEff/snpEff.config GRCh38.99 ../calling_variants/filtering/filtered-vcf.recode.vcf > snpeff.vcf
        	cd ../
	else
        	echo "vcf already annotated"
	fi
}

process geneSelection {

	echo true

	shell:
	'''
	if [ ! -d "annotations/CYP2C19" ]; then
	        echo "selecting CYP2C19 info"
        	mkdir -p annotations/CYP2C19
        	cd annotations/CYP2C19
        	grep CYP2C19 ../snpeff.vcf > CYP2C19.vcf
        	cd ../../
	else
        	echo "CYP2C19 already selected"
	fi
	'''
}

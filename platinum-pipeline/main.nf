#!/home/jovyan/bin/ nextflow

params.trim = false
params.bwaindex = false
params.samtoolsindex = false
params.picarddict = false
params.snpeffgenome = false

channel_name1 = "ch_runfastqc_1"
channel_name2 = "ch_runfastqc_2"

if (params.trim) {
	channel_name1 = "ch_trimlq1"
	channel_name2 = "ch_trimlq2"
}

process getFastq {
	
	echo true

	output:
	file 'SRR1518011_1.fastq.gz' into ch_getfastq_1
	file 'SRR1518011_2.fastq.gz' into ch_getfastq_2

	shell:
	'''
        echo "Downloading SRR1518011_1."
        wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/001/SRR1518011/SRR1518011_1.fastq.gz
        echo "Downloading SRR1518011_2."
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/001/SRR1518011/SRR1518011_2.fastq.gz
	'''
}

process runFastqc {
	
	echo true
	
	input:
	file 'SRR1518011_1.fastq.gz' from ch_getfastq_1
	file 'SRR1518011_2.fastq.gz' from ch_getfastq_2
		
	output:
	file 'SRR1518011_1.fastq.gz' into ch_dummy1
        file 'SRR1518011_2.fastq.gz' into ch_dummy2

	shell:
	'''
	echo "running fastq analysis"
        fastqc SRR1518011_1.fastq.gz SRR1518011_2.fastq.gz
	'''
}

ch_dummy1.view().into { ch_into_trim1; ch_switch1 }
ch_dummy2.view().into { ch_into_trim2; ch_switch2 }


process removePrimersAdapters {
	
	echo true

	when:
    	params.trim

	input:
        file 'SRR1518011_1.fastq.gz' from ch_into_trim1
        file 'SRR1518011_2.fastq.gz' from ch_into_trim2

        output:
        file 'qc/trim/paired/SRR1518011_1.fastq.gz' into ch_removepa_1
        file 'qc/trim/paired/SRR1518011_2.fastq.gz' into ch_removepa_2
	
	shell:
	'''
      	echo "trimming adapters and primers"
      	mkdir -p qc/trim/logs
       	mkdir -p qc/trim/paired
       	mkdir -p qc/trim/unpaired
       	TrimmomaticPE -threads 1 -phred33 -trimlog qc/trim/logs/trimm_logfile \
       	SRR1518011_1.fastq.gz SRR1518011_2.fastq.gz \
       	qc/trim/paired/SRR1518011_1.fastq.gz qc/trim/unpaired/SRR1518011_2.fastq.gz \
       	qc/trim/paired/SRR1518011_2.fastq.gz qc/trim/unpaired/SRR1518011_2.fastq.gz \
       	ILLUMINACLIP:/home/jovyan/coursework-pipeline/resources/primers_adapters.fa:2:30:10 MINLEN:36
	'''
}

process trimLowQuality {

	echo true

	when:
        params.trim

	input:
	file 'SRR1518011_1.fastq.gz' from ch_removepa_1
        file 'SRR1518011_2.fastq.gz' from ch_removepa_2
	
	output:
	file 'qc/trim2/paired/SRR1518011_1.fastq.gz' into ch_trimlq1
	file 'qc/trim2/paired/SRR1518011_2.fastq.gz' into ch_trimlq2
        
	shell:
	'''
	echo "trimming low quality reads"
       	mkdir -p qc/trim2/logs
       	mkdir -p qc/trim2/paired
       	mkdir -p qc/trim2/unpaired
       	TrimmomaticPE -threads 1 -phred33 -trimlog qc/trim2/logs/trimm_logfile \
       	SRR1518011_1.fastq.gz SRR1518011_2.fastq.gz \
       	qc/trim2/paired/SRR1518011_1.fastq.gz qc/trim2/unpaired/SRR1518011_1.fastq.gz \
       	qc/trim2/paired/SRR1518011_2.fastq.gz qc/trim2/unpaired/SRR1518011_2.fastq.gz \
       	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	'''
}

ch_switch1.until{params.trim}.mix(ch_trimlq1).set{ch_into_alignment1}
ch_switch2.until{params.trim}.mix(ch_trimlq2).set{ch_into_alignment2}

process indexReference1 {
	
	echo

        when:
        params.bwaindex

	shell:
	'''
	echo "indexing reference genome"
        bwa index /home/jovyan/coursework-pipeline/resources/reference/reference-genome.fasta
	'''
}

process alignSequence {

	echo true
	
	input:
        file 'SRR1518011_1.fastq.gz' from ch_into_alignment1
        file 'SRR1518011_2.fastq.gz' from ch_into_alignment2
	
	output:
	file 'initial-output.sam' into ch_alignsequence	

	shell:
	'''
        echo "aligning sequence"
        bwa mem -t 8 -R "@RG\\tID:identifier\\tLB:library\\tPL:platform\\tPU:platform-unit\\tSM:sample-name" \
        /home/jovyan/coursework-pipeline/resources/reference/reference-genome.fasta \
       	SRR1518011_1.fastq.gz \
        SRR1518011_2.fastq.gz > initial-output.sam
	'''
}

process convertSam {

	echo true

	input:
	file 'initial-output.sam' from ch_alignsequence

	output:
	file 'initial-output.bam' into ch_convertsam	

	shell:
	'''
        echo "converting sam to bam"
        samtools view -b initial-output.sam -o initial-output.bam
	'''
}

process alignBam {

	echo true

	input:
	file 'initial-output.bam' from ch_convertsam

	output:
	file 'sorted-file.bam' into ch_alignbam
	
	shell:
	'''
        echo "sorting bam"
        samtools sort initial-output.bam -o sorted-file.bam
	'''
}

process indexBam1 {

	echo true

	input:
	file 'sorted-file.bam' from ch_alignbam

	output:
	file 'sorted-file.bam' into ch_indexbam
	file 'sorted-file.bam.bai' into ch_indexbai	
	
	shell:
	'''
        echo "indexing bam"
        samtools index sorted-file.bam
	'''
}

process runFlagstat {

	echo true
	
	input:
	file 'sorted-file.bam' from ch_indexbam
	file 'sorted-file.bam.bai' from ch_indexbai

	output:
	file 'sorted-file.bam' into ch_runflagstat_bam
	file 'sorted-file.bam.bai' into ch_runflagstat_bai

	shell:
	'''
        echo "generating stats"
        samtools flagstat sorted-file.bam > flagstat.txt
	'''
}

process markingDuplicates {

	echo true

	input:
	file 'sorted-file.bam' from ch_runflagstat_bam
	file 'sorted-file.bam.bai' from ch_runflagstat_bai	

	output:
	file 'removed-duplicates.bam' into ch_markingduplicates

	shell:
	'''
        echo "marking duplicates"
        java -jar ~/bin/picard/build/libs/picard.jar MarkDuplicates INPUT=sorted-file.bam OUTPUT=removed-duplicates.bam METRICS_FILE=duplicate-metrics.txt
	'''
}

process indexBam2 {
	
	echo true

	input:
	file 'removed-duplicates.bam' from ch_markingduplicates

	output:
	file 'removed-duplicates.bam' into ch_indexbam2
	file 'removed-duplicates.bam.bai' into ch_indexbai2

	shell:
	'''
        echo "indexing removed duplicates bam file"
        samtools index removed-duplicates.bam
	'''
}

process runQualimap {

	echo true

	input:
        file 'removed-duplicates.bam' from ch_indexbam2
        file 'removed-duplicates.bam.bai' from ch_indexbai2	

	output:
	file 'removed-duplicates.bam' into ch_runqualimap_bam
        file 'removed-duplicates.bam.bai' into ch_runqualimap_bai
	
	shell:
	'''
        echo "producing qualimap qc report"
        qualimap bamqc -bam removed-duplicates.bam -outdir quality-map-report
	'''
}

process indexReference2 {

	echo true

        when:
        params.samtoolsindex

	shell:
	'''
        echo "indexing reference genome"
        samtools faidx /home/jovyan/coursework-pipeline/resources/reference/reference-genome.fasta
	'''
}

process createReferenceDict {

	echo true

        when:
        params.picarddict

	shell:
	'''
        echo "creating reference dictionary"
        java -jar ~/bin/picard/build/libs/picard.jar CreateSequenceDictionary -R /home/jovyan/coursework-pipeline/resources/reference/reference-genome.fasta
	'''
}

process recalibrateData {

	echo true

	input:
	file 'removed-duplicates.bam' from ch_runqualimap_bam
        file 'removed-duplicates.bam.bai' from ch_runqualimap_bai

	output:
	file 'removed-duplicates.bam' into ch_recalibrate_bam
        file 'removed-duplicates.bam.bai' into ch_recalibrate_bai
	file 'recal_data.table' into ch_recalibrate_table

	shell:
	'''
        echo "recalibrating data"
        gatk BaseRecalibrator -I removed-duplicates.bam -R /home/jovyan/coursework-pipeline/resources/reference/reference-genome.fasta \
        --known-sites /home/jovyan/coursework-pipeline/resources/resources-broad-hg38-v0-Homo_sapiens_assembly38.known_indels.vcf \
        --known-sites /home/jovyan/coursework-pipeline/resources/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf -O recal_data.table
	'''
}

process rescoreData {

	echo true

	input:
        file 'removed-duplicates.bam' from ch_recalibrate_bam
        file 'removed-duplicates.bam.bai' from ch_recalibrate_bai
        file 'recal_data.table' from ch_recalibrate_table

	output:
	file 'rescored.bam' into ch_rescored

	shell:
	'''
        echo "apply base quality score recalibration"
        gatk ApplyBQSR -R /home/jovyan/coursework-pipeline/resources/reference/reference-genome.fasta -I removed-duplicates.bam -bqsr recal_data.table -O rescored.bam
	'''
}

process haplotypeCalling {

	echo true

	input:
	file 'rescored.bam' from ch_rescored

	output:
	file 'rescored.bam' into ch_haplotype_bam
	file 'gatk-file.vcf' into ch_haplotype_vcf
	
	shell:
	'''
        echo "haplotype caller running"
        gatk HaplotypeCaller -R /home/jovyan/coursework-pipeline/resources/reference/reference-genome.fasta -I rescored.bam -mbq 20 --minimum-mapping-quality 50 -O gatk-file.vcf
	'''
}

process selectingVariants {
	
	echo true

	input:
	file 'gatk-file.vcf' from ch_haplotype_vcf

        output:
        file 'selected-gatk.vcf' into ch_selectvariants
	
	shell:
	'''
	echo "selecting variants"
	gatk SelectVariants -R /home/jovyan/coursework-pipeline/resources/reference/reference-genome.fasta --variant gatk-file.vcf -O initial-selected-gatk.vcf --select-type SNP
	/home/jovyan/coursework-pipeline/platinum-pipeline/scripts/selection.sh
	'''
}

process runFreebayes {

	echo true

	input:
	file 'rescored.bam' from ch_haplotype_bam

	output:
	file 'freebayes-file.vcf' into ch_freebayes

	shell:
	'''
        echo "running freebayes"
        freebayes -q 20 -m 50 -u -f /home/jovyan/coursework-pipeline/resources/reference/reference-genome.fasta rescored.bam > freebayes-initial-file.vcf
	/home/jovyan/coursework-pipeline/platinum-pipeline/scripts/freebayes.sh
	'''
}

process compareCalling {

	echo true

	input:
	file 'selected-gatk.vcf' from ch_selectvariants
	file 'freebayes-file.vcf' from ch_freebayes

	output:
	file 'freebayes-file.vcf' into ch_compare	

	shell:
	'''
        echo "comparing variant calling"
        vcftools --vcf selected-gatk.vcf --diff freebayes-file.vcf --diff-site --out compare-vcf-files \
	'''
}

process filterVariants {

	echo true

	input:
	file 'freebayes-file.vcf' from ch_compare

	output:
	file 'filtered-vcf.recode.vcf' into ch_filtervariants	

	shell:
	'''
        echo "filter chosen vcf file"
        vcftools --vcf freebayes-file.vcf --minDP 3 --minQ 20 --out temp --recode --recode-INFO-all
        vcftools --vcf temp.recode.vcf --max-missing 1 --out filtered-vcf --recode --recode-INFO-all
	'''
}

process downloadSnpeffGenome {

	echo true

	when:
        params.snpeffgenome

	shell:
	'''
	echo "downloading snpeff genome"
        java -jar ~/bin/snpEff/snpEff.jar download GRCh38.99
       	cp -r ~/bin/snpEff/data/GRCh38.99 /home/jovyan/coursework-pipeline/resources/snpeff-genome/
	'''
}

process annotateVCF {

	echo true

	input:
	file 'filtered-vcf.recode.vcf' from ch_filtervariants

	output:
	file 'snpeff.vcf' into ch_annotate
	
	shell:
	'''
        echo "annotating vcf file"
	java -jar ~/bin/snpEff/snpEff.jar -c ~/bin/snpEff/snpEff.config GRCh38.99 filtered-vcf.recode.vcf > snpeff.vcf
	'''
}

process geneSelection {

	echo true

	input:
	file 'snpeff.vcf' from ch_annotate
	
	shell:
	'''
	echo "selecting CYP2C19 info"
        grep CYP2C19 snpeff.vcf > CYP2C19.vcf
	'''
}

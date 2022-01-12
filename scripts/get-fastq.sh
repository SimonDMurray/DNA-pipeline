#!/bin/bash
set -euo pipefail

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

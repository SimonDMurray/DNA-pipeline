ARG ubuntu_version=18.04

#ENV PATH="/root/miniconda3/bin:${PATH}"
#ARG PATH="/root/miniconda3/bin:${PATH}"

FROM ubuntu:$ubuntu_version

RUN set -xe \
	&& apt-get -y update \
    	&& apt-get -y install fastqc \
    	&& apt-get -y install trimmomatic \
    	&& apt-get -y install libncurses5 \
    	&& apt-get -y install wget \ 
	&& apt-get install -y git \   
	&& apt-get -y install vcftools

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
	&& mkdir /root/.conda \
    	&& bash Miniconda3-latest-Linux-x86_64.sh -b \
    	&& rm -f Miniconda3-latest-Linux-x86_64.sh 

RUN /root/miniconda3/bin/conda install -c bioconda -y samtools

RUN /root/miniconda3/bin/conda install -c bioconda -y bwa

RUN /root/miniconda3/bin/conda install -c bioconda -y freebayes

RUN /root/miniconda3/bin/conda install -c bioconda -y qualimap

RUN /root/miniconda3/bin/conda install -c bioconda -y gatk4

RUN git clone https://github.com/broadinstitute/picard.git \
	&& cd picard \
	&& ./gradlew shadowJar

ENV PATH="/root/miniconda3/bin:${PATH}"

# parent image where image is based on
FROM ubuntu:22.04

# image is non-interactive ONLY DURING BUILD
ARG DEBIAN_FRONTEND=noninteractive

# run updates and install essentials and fix possible broken packages
#RUN apt-get update -y --fix-missing
#RUN apt install -y software-properties-common
#RUN add-apt-repository ppa:deadsnakes/ppa
#RUN apt-get install -y --fix-missing build-essential g++ cmake git-all wget unzip python3.9 curl python3.9-distutils libpython3.9-dev python3.9-venv
#RUN curl https://bootstrap.pypa.io/get-pip.py | python3.9

RUN apt-get update && apt-get install -y --fix-missing build-essential g++ cmake git-all
RUN apt-get update && apt-get install -y \
	wget \
	unzip \
	python3 \
	python3-pip \
	bedtools \
	bowtie2 \
&& rm -rf /var/lib/apt/lists/*

# install primer3 (version 2.6.1)
WORKDIR /bin
RUN wget https://github.com/primer3-org/primer3/archive/v2.6.1.tar.gz
RUN tar -xvf v2.6.1.tar.gz
WORKDIR /bin/primer3-2.6.1/src
RUN make all

WORKDIR /

# install bowtie2
#WORKDIR /bin
#RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.4.5/bowtie2-2.4.5-linux-x86_64.zip
#RUN unzip bowtie2-2.4.5-linux-x86_64.zip
#WORKDIR /

# install bigBedtoBed
RUN apt-get update && apt-get install -y libkrb5-3 \
	curl && rm -rf /var/lib/apt/lists/*
WORKDIR /bin
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
RUN chmod a+x bigBedToBed

WORKDIR /

# install NUPACK
ADD nupack-4.0.1.3 nupack-4.0.1.3
RUN pip3 install -U nupack -f nupack-4.0.1.3/package
RUN pip3 install -U pytest
RUN python3 -m pytest --pyargs nupack

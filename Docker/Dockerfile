# Docker filer for DMAS
# Arne Blom, 2024

# Base image
FROM ubuntu:24.04

# image is non-interactive ONLY DURING BUILD
ARG DEBIAN_FRONTEND=noninteractive

# install 'make' function (in build-essential & others from primer3 manual)
RUN apt-get update -y --fix-missing
RUN apt-get install -y --fix-missing build-essential g++ cmake git-all
RUN apt-get install -y wget
RUN apt-get install -y unzip
RUN apt-get install -y python3
RUN apt-get install -y python3-pip
RUN apt-get install -y git

# install python packages
RUN python3 -m pip config set global.break-system-packages true
RUN python3 -m pip install -U biopython
RUN python3 -m pip install -U bio
RUN python3 -m pip install -U pandas
RUN python3 -m pip install -U datetime
RUN python3 -m pip install -U pybedtools
RUN python3 -m pip install -U requests

# install primer3
WORKDIR /bin
RUN git clone https://github.com/primer3-org/primer3.git primer3
WORKDIR /bin/primer3/src
RUN make all
RUN cp ./primer3_core /usr/bin/primer3_core
WORKDIR /bin
RUN rm -rf primer3

# install bigBedToBed dependencies
RUN apt-get install -y libkrb5-3 curl

# bigBedToBed
RUN wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
RUN chmod +x bigBedToBed

# install bowtie2
RUN apt-get install -y bowtie2


# ViennaRNA (RNAfold)
WORKDIR /bin
RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_6_x/ViennaRNA-2.6.4.tar.gz
RUN tar -xzvf ViennaRNA-2.6.4.tar.gz 
WORKDIR /bin/ViennaRNA-2.6.4
RUN ./configure
RUN make
RUN make install
WORKDIR /bin 
RUN rm -rf ViennaRNA-2.6.4 ViennaRNA-2.6.4.tar.gz

# bedtools
RUN apt-get install -y bedtools
RUN rm -rf /var/lib/a

# scripts
WORKDIR /
ADD scripts /usr/bin
RUN chmod +x /usr/bin/00_Create_input.py /usr/bin/01_generate_template.py /usr/bin/02_check_coord.py /usr/bin/03_get_sec_str_temp_ViennaRNA.py /usr/bin/04_get_SNPs.py /usr/bin/05_primer_generation.py /usr/bin/06_Primer3_common_primer.py /usr/bin/07_primer_validation.py /usr/bin/08_get_sec_str_amp_ViennaRNA.py /usr/bin/09_primer_specificity.py /usr/bin/10_filter.py

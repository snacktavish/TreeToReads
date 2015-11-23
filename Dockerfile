FROM phusion/baseimage:latest
MAINTAINER Justin Payne, justin.payne@fda.hhs.gov

WORKDIR /sw/

RUN apt-get update && apt-get install -y \
		git \
		seq-gen \
		python \
		python-dev \
		python-pip \
		samtools \
	&& curl -O http://www.niehs.nih.gov/research/resources/assets/docs/artbinchocolatecherrycake031915linux64tgz.tgz \
	&& tar -xzvf ./artbinchocolatecherrycake031915linux64tgz.tgz \
	&& pip install dendropy \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*
	
RUN git clone https://github.com/snacktavish/TreeToReads.git

ENV PATH $PATH:/sw/art_bin_ChocolateCherryCake

WORKDIR /sw/TreeToReads/

ENTRYPOINT ["python","/sw/TreeToReads/treetoreads.py"]

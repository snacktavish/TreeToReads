FROM phusion/baseimage:latest
MAINTAINER Justin Payne, justin.payne@fda.hhs.gov

WORKDIR /tmp/
RUN apt-get update && apt-get install -y \
	git \
	seq-gen \
	python \
	python-dev \
	python-pip 
RUN curl -O http://www.niehs.nih.gov/research/resources/assets/docs/artbinvanillaicecream031114linux64tgz.tgz
RUN tar -xzvf ./artbinvanillaicecream031114linux64tgz.tgz
RUN export PATH="$PATH:$PWD/artbinvanillaicecream031114linux64"

RUN pip install dendropy numpy
RUN git clone https://github.com/snacktavish/TreeToReads.git

WORKDIR /tmp/TreeToReads/

ENTRYPOINT python /tmp/TreeToReads/treetoreads.py
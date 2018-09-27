#! /bin/bash

## This script prepares an `isolates.tab` file as an input for a `nullarbor` run on simulated data

# move to TTS output directory
TTS_outdir=$1 
cd ${TTS_outdir}

# prepare isolates.tab
ls -d fastq/* | cut -d'/' -f 2 > isolist
ls fastq/*/*1.fq.gz > read1
ls fastq/*/*2.fq.gz > read2

paste isolist read1 read2 > isolates.tab 

rm isolist read1 read2

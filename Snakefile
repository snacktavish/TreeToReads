"""
A snakefile for simulating a cluster-detection job using TreeToReads, and evaluating with Nullarbor/Snippy
"""

import datetime
import glob
from pathlib import Path
# from Bio import SeqIO
import sys
# import pandas as pd

# genome_file, job_name, job_tree = sys.argv[1], sys.argv[2], sys.argv[3]

# configs variables
genome_file = config['genome']
job_name = config['name']
job_tree = config['tree']

genome_path, tree_path = Path( genome_file ), Path( job_tree )

genome_local = Path( job_name, genome_path.name )
tree_local = Path( job_name, "Reference.fasta" )
config_local = Path( job_name, job_name + ".config" )


def check_genome( Gfile ):
	if Gfile.exists() and Gfile.name.endswith(".fasta") or Gfile.name.endswith(".fa"):
		return( Gfile )
	elif Gfile.exists() and Gfile.name.endswith(".gbk"):
   		myFasta = Gfile.name + ".fasta"
		SeqIO.convert( Gfile, "genbank", myFasta, "fasta" )
		return( myFasta )


# rule all:
#     input: 
#         report=all_target(config['mdu_report']),
#         input_file=config['input_file']

rule get_genome:
	input:
		genome_file
	output:
		str( genome_local )
	shell:
		"cp {input} {output}"

rule get_newick:
    input:
        job_tree
    output:
        str( tree_local )
    shell:
        "cp {input} {output}"

rule write_config:
    input:
        str( tree_local ),
        str( genome_local )
    output:
        str( config_local )
    shell:
        "python write_config.py {input[1]} %s {input[0]}" % job_name

rule run_TTR:
    input:
        config_local
    output:
        job_name + "/fastq/sim_Reference/sim_Reference_1.fq.gz",
        job_name + "/fastq/sim_Reference/sim_Reference_2.fq.gz"
    shell:
        "python treetoreads.py {input}"




#! /bin/python3

## This script is intended to accompany my Snakefile for TreeToReads -- it should take inline arguments and write them into a *.config file

import sys

genome_file = sys.argv[1]
job_name = sys.argv[2]
job_tree = sys.argv[3]

#print( genome_file, job_name, job_tree )

N_var_sites = 10
coverage = 60


REQUIRED_PARAMS = f""" 
treefile_path = {job_tree}
number_of_variable_sites = {N_var_sites}
base_genome_name = Reference
base_genome_path = {genome_file}
output_dir = {job_name}
prefix = sim_
"""

EVO_PARAMS = f"""
rate_matrix = 1,1,1,1,1,1
mutation_clustering = ON
percent_clustered = 0.25
exponential_mean = 25
gamma_shape = 5
"""

OPTION_PARAMS = f"""
coverage = {coverage}
read_length = 150
fragment_size = 380
stdev_frag_size = 120
"""

#print( REQUIRED_PARAMS, EVO_PARAMS, OPTION_PARAMS )

# ART Optional parameters (for more fine grained control ART can be run seperately on the mutated genomes found in outdir/fasta_files)
# error_model1 = example/ErrprofR1.txt  # If you haven't generated have one of your own using ART, you can use one supplied by ART.
# error_model2 = example/ErrprofR2.txt  # Un-comment these lines (delete the first #) to set a non-default error profile

with open( "%s/%s.config" % ( job_name, job_name), 'w' ) as outfile:
	outfile.write( REQUIRED_PARAMS + '\n' )
	outfile.write( EVO_PARAMS + '\n' )
	outfile.write( OPTION_PARAMS )
	outfile.close()

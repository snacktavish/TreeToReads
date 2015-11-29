Tree to reads tutorial
======================

Tree to reads (TTR) is a simulation pipeline to generate next generation sequencing reads from realistic phylogenies.
Can be used to test effects of model of evolution, rates of evolution, genomic distribution of mutations, and phylogenetic relatedness of samples and of reference genome on SNP calling and evolutionary inference.

Inputs are a phylogeny, a genome to be used as a tip in the tree, and a set of configuration parameters in a control file

Optional inputs can include a sequencing error model parameterized from empirical data, and a distribution for the distances separating pairs of mutations in the genome.

Outputs are mutated genomes representing all tips in the phylogeny, and simulated whole genome sequencing reads representing those genomes. These are are useful for testing and comparison of analysis pipelines. Mutations are currently only single nucleotide variants - no indels or rearrangements.

The code is still in draft format - but testing welcome, and will be supported via email ejmctavish, gmail. 


## Configuration
All the necessary parameters for TTR are specified in a configuration file. In this tutorial we will walk though the parameter arguments in the "ttr.cfg" example file.

* Choose a phylogeny (we will use example/simtree.tre). Include the full path to the tree file in either newick or nexus in the configuration file.  
```
treefile_path = example/simtree.tre
```
NOTE: The branch lengths in this tree will be proportional to the branch lengths in your outputs, but not equivalent, as the final branch lengths will depend on the number of variable sites selected. Long branch lengths in the input tree will result in more multiple hit mutations at the same sites.  
ALSO: TTR cannot handle polytomies, so any polytomies in this tree will be randomly resolved with 0 length branch lengths using dendropy, and saved to the output dir as simtree.tre

* Select a "base genome" and specify the path to this file in the configuration file. This the base genome on which mutations will be placed. Currently TTR will only work with files containing a single contig or chromosome. Simulate reads for separate chromosomes individually.
One tip in your final tree will have this genome sequence. This tip label should be specified using the parameter "base_genome_name".  
The number of variable sites can vary between 0 and the length of the genome.
```
base_genome_name = gi  
base_genome_path = example/mini_ref.fasta  
number_of_variable_sites = 20
```
* Specify an output directory to which files should be written. Take care, because output will overwrite previous runs in that directory
```
output_dir = example_out
```
* Set the parameters of your evolutionary model. These parameters may be generated from real sequence analysis, or a null model. TTR requires a full general time reversible evolutionary model including a rate matrix, state frequencies, and a gamma shape parameter. These parameters are passed to seq-gen (http://bioweb2.pasteur.fr/docs/seq-gen/) and are in seq-gen like format.
A good ways to estimate these model parameters from your observed SNP data is using Paup (http://paup.csit.fsu.edu/) and ModelTest (http://www.molecularevolution.org/software/phylogenetics/modeltest).
```
rate_matrix = 1,1,1,1,1,1
```
where these 6 values are decimal numbers for the relative rates of substitutions from (for nucleotides) A to C, A to G, A to T, C to G, C to T and G to T respectively, separated by commas.
```
freq_matrix =  0.19,0.31,0.29,0.22
```
where these 4 values are decimal numbers for the frequencies of the nucleotides A, C, G and T
```
gamma_shape = 5
```
where the value is a real number >0 that specifies the shape of the gamma distribution to use with gamma rate heterogeneity.

* Set the parameters for read simulation.  
TTR currently only automates simulation of illumina paired end short read data (but that will be expanded in future).  The read simulation parameters are passed to Art (http://www.niehs.nih.gov/research/resources/software/biostatistics/art/).
For more fine grained control, or simulation of single end illumina 454 or SOLiD reads, ART can be run separately on the mutated genomes found in outdir/fasta_files, using the instructions found in the ART/README.  
The only required parameter is average coverage across the genome.
```
coverage = 20
```
To run TTR to generate simulated genomes at the tips of the tree, but not generate simulated reads from these genomes, set coverage to none.
```
coverage = 0
```
* Optional read simulation parameters:  
You will likely want to change the read lengths and fragment sizes.
The defaults are:
```
read_length = 150
fragment_size = 350
stdev_frag_size = 130
```
the maximum possible value for read length is currently 250 bp, under the default error models, but you can build an error model for different read lengths as described below.
ART by default selects a built-in quality score profile according to the read length specified for the run.
To use a default profile 
or you can provide an error model for the reads. You can either build one from your own data using Art or use one of the error models packaged with Art, found in the ART/Illumina_profiles directory. To read more about the Art error models see art_illumina_README in the Art directory.
To build an error model based on your data use:  
```
./Illumina_readprofile_art out_profile_name input_fastq_dir
```
where out_profile_name is the name for your error model, and input_fastq_dir is a directory of observed illumina sequence reads.  
The example uses an error model built from observed data.
```
(in the config file)
error_model1 = example/ErrprofR1.txt  
error_model2 = example/ErrprofR2.txt
```

* Optional mutation distributions.  
By default the location of variable sites will be selected from a random uniform distribution across the genome. Optionally, some proportion of mutation locations can be generated in a clustered fashion (pairs of mutations closer to each other than would be expected from a random distributions). The distance between these sites is drawn from an exponential distribution, the variance of which must be given using exponential_mean. The lower the mean value the more closely clustered paired sites will be. 
These can be tricky to estimate from observed data, and trial and error comparing observed and simulated mutation distance distributions is helpful, but note that SNP calling 
```
#parameters for clustering of variable site locations (OPTIONAL)
mutation_clustering = ON
percent_clustered = 0.25
exponential_mean = 10
```
Setting mutation_clustering clustering to OFF will cause these values to be ignored and mutation locations to be drawn from a uniform random distribution.

## Run the program!
```
python treetoreads.py
```
TTR will automatically look for and read a configuration file named "ttr.cfg".
You can also save your configuration under a different name and pass it as an argument, e.g.
```
python treetoreads.py ttr_alt.cfg
```

## Output files
By default the sequence files will have the names of the tips in the input tree.
Alternatively, a prefix can be specified using
```
prefix = sim_
```
The key output files consist of:
    * fasta_files - a folder containing the simulated genomes for each tip in the tree  The fasta files in this folder can be used in conjunction with Art for more fine grained read simulation.  
    * fastq - folder containing folders with the names of each tip from the simulation tree. In each of these folders are the simulated reads  in .fastq format and .sam format file of the read alignments.  
    * mutsites.txt - unordered list of the locations of mutations in the genome  
    * var_site_matrix - an unordered list of the base present in each tip at each variable site, in the format "tip_name base genome_location"  

Other files
    * analysis_configuration.cfg is a copy of the config file used for the analysis  
    * seq_sim.txt is the seqgen output file from which the variable sites are drawn.  
    * art_log and seqgen_log are the output files of art and seqgen, respectively, and are useful for diagnosing issues.  
    * simtree.tre and simtree.tre.bu are the input tree with any polymotimies resolved.  

#TreeToReads

Simulation pipeline to generate next generation sequencing reads from realistic phylogenies.  
Can be used to test effects of model of evolution, rates of evolution, 
genomic distribution of mutations, and phylogenetic relatedness of samples and of reference genome 
on SNP calling and evolutionary inference.  

Still in draft format - but testing welcome, and will be supported via email ejmctavish, gmail.  

##Requirements:

Installed globally
-   Seq-Gen
-   Art
(can be run without Art if you want to generate mutated genomes, but not reads)

python packages
-   Dendropy


-------------------------

##To install requirements

###Install seq-gen, software to simulate mutations (http://tree.bio.ed.ac.uk/software/seqgen/)
    sudo apt-get install seq-gen

###Install ART

    wget http://www.niehs.nih.gov/research/resources/assets/docs/artbinvanillaicecream031114linux64tgz.tgz
    tar -xzvf artbinvanillaicecream031114linux64tgz.tgz

You need to add art_illumina to your path


###Install Dendropy

    easy_install dendropy


-----------------------------------------------------------
##Running the simulations:

    git clone https://github.com/snacktavish/TreeToReads.git
    cd TreeToReads
    python treetoreads.py seqsim.cfg
 

Edit config file, seqsim.cfg, to fit your data.
The script by default look for a file called 'seqsim.cfg'
or first argument can be the path to a control file.

Currently only runs art_illumina and generates paired end illumina data.
Alternatively, genomes can be generated, and ART run seperately using any chosen parameters.

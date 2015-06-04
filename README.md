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
(can be run with out Art if you want to generate mutated genomes, but not reads)

python packages
-   Dendropy
-   Numpy


-------------------------

##To install requirements


    sudo apt-get install seq-gen
    wget http://www.niehs.nih.gov/research/resources/assets/docs/artbinvanillaicecream031114linux64tgz.tgz
    tar -xzvf artbinvanillaicecream031114linux64tgz.tgz
So, then you need to add artbinvanillaicecream031114linux64tgz to your path. A bit annoying.

    easy_install dendropy

-----------------------------------------------------------
##Running the simulations:

    git clone https://github.com/snacktavish/TreeToReads.git
    cd TreeToReads

Edit config file, seqsim.cfg, to fit your data!

    python simulations.py mysims.cfg

by default looks for seqsim.cfg, or first argument can be a control file

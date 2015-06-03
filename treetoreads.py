#!/usr/bin/python
import dendropy
from subprocess import call
import random
import os
import sys
import argparse
import numpy
import sys

VERSION="0.0.1"

class TreeToReads:
  _argread=0
  _treeread=0
  _simran=0
  _madeout=0
  _siteread=0
  _mutlocs=0
  _genmut=0
  _genread=0
  _vargen=0
  def __init__(self, configfi='seqsim.cfg', run=1):
        """A method to read a tree, resolve polytomes, generate mutations and simulate reads."""
        self.configfi=configfi
        if os.path.isfile(self.configfi):
          sys.stdout.write("Running TreetoReads using configuration file {}\n".format(self.configfi))
        else:
          sys.stderr.write("Config file '{}' not found. Exiting.\n".format(self.configfi))
          sys.exit()
        self.run=run
        if self.run:
          self.runSims()
        self._checkArgs()
#        test_deps()
  def test_deps(self):
    t=1
  def runSims(self):
            self.runART()      
  def readArgs(self):
    self._argread=1
    if len(sys.argv) > 1:
      self.configfi=sys.argv[1]
      try:
        config=open(self.configfi).readlines()
      except:
        sys.stderr.write("Couldn't open config file {}. Exiting.\n".format(self.configfi))
        sys.exit()
    else:
      try:
        config=open(self.configfi)
      except:
        sys.stderr.write("Config file '{}' not found. Exiting.".format(self.configfi))
        sys.exit()
    self.config={}
    for lin in config:
      lii=lin.split('=')
      self.config[lii[0].strip()]=lii[-1].split('#')[0].strip()
    self._checkArgs()
    sys.stdout.write("Arguments read\n")
  def makeOut(self):
      if not self._argread:
        self.readArgs()
      self._madeout=1
      sys.stdout.write('output directory is {}\n'.format(self.outd))
      if not os.path.isdir(self.outd):
          os.mkdir(self.outd)
      self.bashout = open('{}/analysis.sh'.format(self.outd),'w')
      configout = open('{}/analysis_configuration.cfg'.format(self.outd),'w')
      for lin in open(self.configfi).readlines():
          configout.write('lin')
      configout.close()
  def getArg(self,nam):
    try:
      len(self.argdict)
    except:
      self.readArgs()
    if nam in self.argdict:
        return(self.config[self.argdict[nam]])
    else:
      try:
        return(self.config[nam])
      except:
        return(None)
  def _checkArgs(self):
    if self._argread!=1:
      self.readArgs()
    self.argdict={'treepath':'treefile_path', 'nsnp':'number_of_snps', 'anchor_name':'anchor_name', 'genome':'anchor_genome_path', 'ratmat':'rate_matrix', 'freqmat':'freq_matrix', 'shape':'shape', 'errmod1':'error_model1', 'errmod2':'error_model2', 'cov':'coverage', 'outd':'output_dir'}
    for arg in self.argdict:
      if self.argdict[arg] not in self.config:
        sys.err.write("{} is missing from the config file".format(self.argdict[arg]))
        sys.exit()
    try:
        self.nsnp=int(self.getArg('nsnp'))
        sys.stdout.write('Number of SNPS is {}\n'.format(self.nsnp))
    except:
        sys.stderr.write("number of SNPs {} could not be coerced to an integer. Exiting.\n".format(self.getArg('nsnp')))
        sys.exit()
    if not len(self.getArg('ratmat').split(','))==6:
        sys.stderr.write("{} values in rate matrix, there should be 6. Exiting.\n".format(len(self.getArg('ratmat').split(','))))
        sys.exit()
    if not len(self.getArg('freqmat').split(','))==4:
        sys.stderr.write("{} values in freq matrix, there should be 4. Exiting.\n".format(len(self.getArg('freqmat').split(','))))
        sys.exit()  
    try:
        float(self.getArg('shape'))
    except:
        sys.stderr.write("shape parameter {} could not be coerced to a float. Exiting.\n".format(self.getArg('shape')))
        sys.exit() 
    try:
        open(self.getArg('treepath'))
    except:
      sys.stderr.write("cCuld not open treefile {}. Exiting.\n".format(self.getArg('treepath')))
      sys.exit() 
    try:
        open(self.getArg('genome'))
    except:
        sys.stderr.write("Could not open anchor genome {}. Exiting.\n".format(self.getArg('genome')))
        sys.exit() 
    try:
        self.outd=self.getArg('outd')
    except:
        self.outd=('ttr_out')
    if self.getArg('mutation_clustering') is not None:
      if self.getArg('mutation_clustering')=='ON':
        self.clustering=1
        try:
            self.clustPerc = float(self.getArg('percent_clustered'))
            self.lambd = float(self.getArg('exponential_lambda'))
            print("clustering proprotion is {}".format(self.clustPerc))
            print("lambda is {}".format(self.lambd))
        except:
            sys.sterr.write("Problem reading clustering parameters, requires number for 'percent_clustered' and 'exponential_lambda. Exiting.'\n")
            sys.exit()
      else: 
        sys.stdout.write('Mutation clustering is OFF, to use set mutation_clustering = ON and values for "percent_clustered" and "exponential_lambda"\n')
        self.clustering=0
    else: 
        sys.stdout.write('Mutation clustering is OFF\n')
        self.clustering=0
  def readTree(self):
    self._treeread=1
    if not self._madeout:
      self.makeOut()
    #import tree from path
    taxa=dendropy.TaxonSet()
    tree=dendropy.Tree.get_from_path(self.getArg('treepath'),'newick',taxon_set=taxa)
    self.seqnames = taxa.labels()
    if not self.getArg('anchor_name') in self.seqnames:
      sys.sterr.write("anchor genome name {} is not in tree. Exiting.\n".format(self.getArg('anchor_name')))
      sys.exit()
    tree.resolve_polytomies()
    if tree.length >= 1:
      sys.stderr.write("WARNING: Tree length is high - scale down tree or expect high multiple hits/homoplasy!\n")
    self.outtree="{}/simtree.tre".format(self.outd)
    tree.write(open(self.outtree,'w'),'newick',suppress_internal_node_labels=True)
    linrun="sed -i.bu -e's/\[&U\]//' {}".format(self.outtree)
    self.bashout.write(linrun+'\n')
    os.system(linrun)
    treefi=open(self.outtree).readlines()
    sys.stdout.write("Tree read\n")
  def readGenome(self):
    self._genread=1
    if not self._argread: self.readArgs()
    genfas=open(self.getArg('genome')).readlines()
    crop=[lin[:-1] for lin in genfas[1:]]
    self.gen="".join(crop)
    self.genlen=len(self.gen)
    sys.stdout.write("Genome has {} bases\n".format(self.genlen))
  def generateVarsites(self):
    self._vargen=1
    ## TODO make model variable, currently only GTR is possible
    if not self._treeread:
      self.readTree()
    self.simloc="{}/seqs_sim.txt".format(self.outd)
    lenseqgen=40*self.nsnp ##TODO - scale to branch lengths?
    seqgenpar=['seq-gen', '-l{}'.format(lenseqgen), '-n1', '-mGTR', '-a{}'.format(self.getArg('shape')), '-r{}'.format(self.getArg('ratmat')), '-f{}'.format(self.getArg('freqmat')), '-or']
    seqcall = call(seqgenpar, stdout=open('{}'.format(self.simloc),'w'), stderr=open('{}/seqgen.out'.format(self.outd),'w'), stdin=open('{}'.format(self.outtree)))
    self.bashout.write(" ".join(seqgenpar + ['<', '{}'.format(self.outtree),'>', '{}/{}'.format(self.outd,self.simloc), '2>', '{}/seqgen.out'.format(self.outd)])+'\n')
    sys.stdout.write("Variable sites generated using seq-gen\n")
  def readVarsites(self):
    self._siteread=1
    if not self._vargen:
      self.generateVarsites()
    nucsets={} 
    with open(self.simloc) as f:
        next(f)
        for lin in f:
            bases=lin.split()[1].strip()
            for i, nuc in enumerate(bases):
              if i not in nucsets:
                nucsets[i]=set()
              nucsets[i].add(nuc)
    varsites=[]
    tb=0
    vs=0
    for i in nucsets:
      if len(nucsets[i]) >= 2: 
        varsites.append(i)
        vs+=1
      if len(nucsets[i]) > 2: 
        tb+=1
    sys.stderr.write('WARNING: {} SNP sites with more than 2 bases out of {} sites. Scale down tree length if this is too high.\n'.format(tb,vs))
    simseqs={}
    with open(self.simloc) as f:
        next(f)
        for lin in f:
          seq=lin.split()[0]
          simseqs[seq]=[]
          bases=lin.split()[1].strip()
          for i in varsites:
                simseqs[seq].append(bases[i])
    try:
      assert(set(self.seqnames)==set(simseqs.keys()))
    except:
      sys.sterr.write(self.seqnames)
      sys.sterr.write(simseqs.keys())
      sys.exit()
    ref=simseqs[self.getArg('anchor_name')]
    self.sitepatts={}
    for nuc in ['A','G','T','C']:
       self.sitepatts[nuc]=[]
    for i, nuc in enumerate(ref):
      site={}
      nucs=set()
      for srr in simseqs:
          site[srr]=simseqs[srr][i]
          nucs.add(simseqs[srr][i])
      assert(len(nucs)>1)
      self.sitepatts[nuc].append(site)  #PICKLE THIS SOMEHOW?!?!
    sys.stdout.write("Variable sites read\n")
  def selectMutsites(self):
    if not self._madeout: 
      self.makeOut()
    if not self._genread: 
      self.readGenome()
    self.mutsite="{}/mutsites.txt".format(self.outd)
    fi=open(self.mutsite,"w")
    rands=set()
    if self.clustering:
        nclust=int((self.nsnp*self.clustPerc)/2)
        ranpairA=random.sample(range(self.genlen),nclust)
        for site in ranpairA:
            rands.add(site)
            diff = 0 #1+exponential ... 
            while diff==0: #RISKY at HIGH LAMBDA!
                diff=int(random.expovariate(self.lambd))
            if (random.choice([0, 1]) or (site-diff < 0)) and (site+diff < self.genlen):
              ranpairB=site+diff
            else:
              ranpairB=site-diff
            rands.add(ranpairB)
        ransingle=random.sample(range(self.genlen),int(self.nsnp*(1-self.clustPerc)))
        rands.update(ransingle)
        for site in rands:
                fi.write(str(site)+'\n')
    else:
            ran=random.sample(range(self.genlen),self.nsnp) # Makes sure this works..., test if is sampling with replacement or no. Prefer W/O replacement.
            rands=set(ran)
            for site in rands:
                fi.write(str(site)+'\n')
    sys.stdout.write("realized number of mutations is {}\n".format(len(rands)))
    self.mutlocs=rands
    # set of random numbers determining where the SNPs really fall. can be drawn from SNPlist or just random.
    fi.close()
    self._mutlocs=1
  def mutGenomes(self):
#    if not self._siteread: self.readVarsites()
    if not self._mutlocs: 
      self.selectMutsites()
    if not self._siteread:
      self.readVarsites()
    self.mut_genos={}
    matout=open("{}/SNPmatrix".format(self.outd),'w')
    patnuc={}
    ri=0
    patnuc['A']=0
    patnuc['G']=0
    patnuc['T']=0
    patnuc['C']=0
    snpdic={}
    for nuc in self.gen:
                ri+=1
                if ri in self.mutlocs:
                  patnuc[nuc]+=1
                  snpdic[ri]=patnuc[nuc]
    for seq in self.seqnames:
        self.mut_genos[seq]=[]
        sys.stdout.write("writng genome for {}\n".format(seq))
        genout=open("{}/sim_{}.fasta".format(self.outd,seq),'w')
        ii=0
        genout.write(">SIM_{}".format(seq))
        for nuc in self.gen:
                if ii%70==0:
                   genout.write('\n')
                ii+=1
                if ii in self.mutlocs:             
                   patt=self.sitepatts[nuc][snpdic[ii]]
                   genout.write(patt[seq])
                   self.mut_genos[seq].append(patt[seq])
                   matout.write("{} {} {}\n".format(seq, patt[seq], ii))
                else:
                    genout.write(nuc)
        genout.write('\n')
        genout.write('\n')
        genout.close()
    matout.close()
    self._genmut=1
    sys.stdout.write("Mutated genomes\n")
  def runART(self):
      if not self._genmut:
            self.mutGenomes()
      for seq in self.seqnames:
        artparam=['art_illumina', '-1', self.getArg('errmod2'), '-2', self.getArg('errmod2'), '-p', '-sam', '-i', '{}/sim_{}.fasta'.format(self.outd,seq), '-l', '150', 
                      '-f', self.getArg('cov'), '-m', '350', '-s', '130', '-o', '{}/sim_{}_'.format(self.getArg('outd'),seq)]
        self.bashout.write(' '.join(artparam)+'\n')
        call(artparam)   
      sys.stdout.write("ART generated reads\n")
      """
      #'-l', '150', <- average read length
      # '-f', '20', <-need to get fold coverage
      #'-m', '350', <- fragment size?? need to check
      #'-s', '130', <- stdev fragment size...
      # pull a random number fromt the length of the genome,
      #art_illumina -1 fullprof.txt -2 fullprof.txt -p -sam -i SIM_.fasta -l 150 -f 20 -m 2050 -s 50 -o sim_test
      #Map those mutations back onto the genomes....
      #Get distributaion of SNP locations from real data?!!
      # simulate reads for each of those genomes.
      """


parser=argparse.ArgumentParser(
    description='''Tree to Reads - A python script to to read a tree, resolve polytomes, generate mutations and simulate reads.''',
    epilog="""Still in development - email ejmctavish@gmail.com with questions, suggestions, issues etc.""")
parser.add_argument("configfi",  nargs='?', default="seqsim.cfg", type=str, help="configuration file path. Optional, defaults to seqsim.cfg")
parser.add_argument('-V', '--version', 
                    action='version',                    
                    version='Tree to reads version {}'.format(VERSION))
args = parser.parse_args()



if __name__ == "__main__":
  ttr=TreeToReads(configfi=args.configfi)

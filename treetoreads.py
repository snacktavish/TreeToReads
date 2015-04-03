#!/usr/bin/python
import dendropy
from subprocess import call
import random
import os
import sys
import argparse
import numpy
import sys


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
        self.run=run
        if self.run:
          self.runSims()
        self._checkArgs()
  def runSims(self):
            self.runART()      
  def readArgs(self):
    print("reading arguments")
    self._argread=1
    if len(sys.argv) > 1:
      self.configfi=sys.argv[1]
      try:
        config=open(self.configfi).readlines()
      except:
        print("couldn't open config file {}".conf)
        sys.exit()
    else:
      try:
        config=open(self.configfi)
      except:
        print("config file '{}' not found".format(self.configfi))
        sys.exit()
    self.config={}
    for lin in config:
      lii=lin.split('=')
      self.config[lii[0].strip()]=lii[-1].split('#')[0].strip()
    self._checkArgs()
  def makeOut(self):
      if not self._argread:
        self.readArgs()
      self._madeout=1
      print('output directory is {}'.format(self.outd))
      if not os.path.isdir(self.getArg('outd')):
          os.mkdir(self.getArg('outd'))
      self.bashout = open('{}/analysis.sh'.format(self.getArg('outd')),'w')
      configout = open('{}/analysis_configuration.cfg'.format(self.getArg('outd')),'w')
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
    print("checking args")
    if self._argread!=1:
      self.readArgs()
    self.argdict={'treepath':'treefile_path', 'nsnp':'number_of_snps', 'anchor_name':'anchor_name', 'genome':'anchor_genome_path', 'ratmat':'rate_matrix', 'freqmat':'freq_matrix', 'shape':'shape', 'errmod1':'error_model1', 'errmod2':'error_model2', 'cov':'coverage', 'outd':'output_dir'}
    for arg in self.argdict:
      if self.argdict[arg] not in self.config:
        print("{} is missing from the config file".format(self.argdict[arg]))
        sys.exit()
    try:
        self.nsnp=int(self.getArg('nsnp'))
        print('Number of SNPS is {}'.format(self.nsnp))
    except:
        print("number of SNPs {} could not be coerced to an integer".format(self.getArg('nsnp')))
        sys.exit()
    if not len(self.getArg('ratmat').split(','))==6:
        print("{} values in rate matrix, there should be 6".format(len(self.getArg('ratmat').split(','))))
        sys.exit()
    if not len(self.getArg('freqmat').split(','))==4:
        print("{} values in freq matrix, there should be 4".format(len(self.getArg('freqmat').split(','))))
        sys.exit()  
    try:
        float(self.getArg('shape'))
    except:
        print("shape parameter {} could not be coerced to a float".format(self.getArg('shape')))
        sys.exit() 
    try:
        open(self.getArg('treepath'))
    except:
      print("could not open treefile {}".format(self.getArg('treepath')))
      sys.exit() 
    try:
        open(self.getArg('genome'))
    except:
        print("could not open anchor genome {}".format(self.getArg('genome')))
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
        except:
            print("Problem reading clustering parameters, requires number for 'percent_clustered' and 'exponential_lambda'")
            sys.exit()
      else: 
        print('Mutation clustering is OFF, to use set mutation_clustering = ON and values for "percent_clustered" and "exponential_lambda"')
        self.clustering=0
    else: 
        print('Mutation clustering is OFF')
        self.clustering=0
  def readTree(self):
    print("read tree")
    self._treeread=1
    if not self._madeout:
      self.makeOut()
    #import tree from path
    taxa=dendropy.TaxonSet()
    tree=dendropy.Tree.get_from_path(self.getArg('treepath'),'newick',taxon_set=taxa)
    self.seqnames = taxa.labels()
    if not self.getArg('anchor_name') in self.seqnames:
      print("anchor genome name {} is not in tree!!".format(self.getArg('anchor_name')))
      sys.exit()
    tree.resolve_polytomies()
    self.outtree="{}/simtree.tre".format(self.getArg('outd'))
    tree.write(open(self.outtree,'w'),'newick',suppress_internal_node_labels=True)
    linrun="sed -i.bu -e's/\[&U\]//' {}".format(self.outtree)
    self.bashout.write(linrun+'\n')
    os.system(linrun)
  def readGenome(self):
    self._genread=1
    if not self._argread: self.readArgs()
    genfas=open(self.getArg('genome')).readlines()
    crop=[lin[:-1] for lin in genfas[1:]]
    self.gen="".join(crop)
    self.genlen=len(self.gen)
    print("Genome has {} bases".format(self.genlen))
  def generateVarsites(self):
    print("gen varsites")
    self._vargen=1
    ## TODO make model variable
    if not self._treeread:
      self.readTree()
    self.simloc="{}/seqs_sim.txt".format(self.getArg('outd'))
    lenseqgen=10*int(self.getArg('nsnp'))
    seqcall=" ".join(['seq-gen', '-l{}'.format(lenseqgen), '-n1', '-mGTR', '-a{}'.format(self.getArg('shape')), '-r{}'.format(self.getArg('ratmat')), '-f{}'.format(self.getArg('freqmat')), '-or','<', '{}'.format(self.outtree),'>', '{}'.format(self.simloc)])
    os.system(seqcall)
    self.bashout.write(seqcall +'\n')
  def readVarsites(self):
    print("read varsite")
    self._siteread=1
    if not self._vargen:
      self.generateVarsites()
    nucsets={} 
    with open(self.simloc) as f:
        next(f)
        for lin in f:
            bases=lin.split()[1].strip()
       #     print("seq {} has bases {}".format(lin.split()[0],bases))
            for i, nuc in enumerate(bases):
              if i not in nucsets:
                nucsets[i]=set()
              nucsets[i].add(nuc)
    varsites=[]
    for i in nucsets:
      if len(nucsets[i]) == 2: #THIS LIMITS TO NO MULT HITS!! DO I want that?
        varsites.append(i)
   # print(varsites)
    simseqs={}
    with open(self.simloc) as f:
        next(f)
        for lin in f:
          seq=lin.split()[0]
 #         print("seq name is {}".format(seq))
          simseqs[seq]=[]
          bases=lin.split()[1].strip()
          for i in varsites:
                simseqs[seq].append(bases[i])
    try:
      assert(set(self.seqnames)==set(simseqs.keys()))
    except:
      print(self.seqnames)
      print(simseqs.keys())
      sys.exit()
    ref=simseqs[self.getArg('anchor_name')]
    print(len(ref))
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
  def selectMutsites(self):
    print("select mutsites")
    if not self._madeout: 
      self.makeOut()
    if not self._genread: 
      self.readGenome()
    self.mutsite="{}/mutsites.txt".format(self.getArg('outd'))
    fi=open(self.mutsite,"w")
    rands=set()
    if self.clustering:
        ranpairA=random.sample(range(self.genlen),(self.nsnp*self.clustPerc)/2)
        rands.add(set(ranpairA))
        for site in ranpairA:
            diff = 0 #1+exponential ... 
            while diff==0: #RISKY at HIGH LAMBDA!
                diff=int(random.expovariate(self.lambd))
            if (random.choice([0, 1]) or (site-diff < 0)) and (site+diff < self.genlen):
              ranpairB=site+diff
            else:
              ranpairB=site-diff
            rands.add(ranpairB)
        ransingle=random.sample(range(self.genlen),self.nsnp*(1-self.clustPerc))
        rands.add(set(ransingle))# TODO check this is legal
        for site in rands:
                fi.write(str(site)+'\n')
    else:
            ran=random.sample(range(self.genlen),self.nsnp) # Makes sure this works..., test if is sampling with replacement or no. Prefer W/O replacement.
            rands=set(ran)
            for site in rands:
                fi.write(str(site)+'\n')
    print("realized number of mutations is {}".format(len(rands)))
    self.mutlocs=rands
    # set of random numbers determining where the SNPs really fall. can be drawn from SNPlist or just random.
    fi.close()
    self._mutlocs=1
  def mutGenomes(self):
    print("run mut genomes")
#    if not self._siteread: self.readVarsites()
    if not self._mutlocs: 
      self.selectMutsites()
    if not self._siteread:
      self.readVarsites()
    self.mut_genos={}
    matout=open("{}/SNPmatrix".format(self.getArg('outd')),'w')
    for seq in self.seqnames:
        self.mut_genos[seq]=[]
        print(seq)
        genout=open("{}/sim_{}.fasta".format(self.getArg('outd'),seq),'w')
        patnuc={}
        patnuc['A']=0
        patnuc['G']=0
        patnuc['T']=0
        patnuc['C']=0
        ii=0
        genout.write(">SIM_{}\n".format(seq))
#        print("current seq is {}".format(seq))
        for nuc in self.gen:
                if ii%70==0:
                   genout.write('\n')
                ii+=1
                if ii in self.mutlocs:
#                   print("ref nuc is {}".format(nuc))
                   patt=random.choice(self.sitepatts[nuc]) # risky at nigh num taxa? or too few varsites - all will have same pattern... But won't run out of sites...
#                   print("selected patt with correct ref nuc is")
#                   print(patt)
                   genout.write(patt[seq])
#                   print("and the base for seq {} with reference base {} at site {} is {}\n".format(seq,nuc,ii,patt[seq]))
                   self.mut_genos[seq].append(patt[seq])
                   matout.write("{} {} {}\n".format(seq, patt[seq], ii))
                else:
                    genout.write(nuc)
        genout.write('\n')
        genout.close()
    matout.close()
    self._genmut=1
  def runART(self):
      print("runART")
      if not self._genmut:
            self.mutGenomes()
       #   print(' '.join(['art_illumina', '-1', '../simB/fullprofR1.txt', '-2', '../simB/fullprofR2.txt', '-p', '-sam', '-i', 'sim_{}.fasta'.format(seq), '-l', '150', '-f', '20', '-m', '350', '-s', '130', '-o', 'sim_{}_'.format(seq)]))
       #   call(['art_illumina', '-1', args['error_model1'], '-2', args['error_model2'], '-p', '-sam', '-i', 'sim_{}.fasta'.format(seq), '-l', '150', '-f', args['coverage'], '-m', '350', '-s', '130', '-o', 'sim_{}_'.format(seq)])
      for seq in self.seqnames:
        artparam=' '.join(['art_illumina', '-1', self.getArg('errmod2'), '-2', self.getArg('errmod2'), '-p', '-sam', '-i', '{}/sim_{}.fasta'.format(self.getArg('outd'),seq), '-l', '150', 
                      '-f', self.getArg('cov'), '-m', '350', '-s', '130', '-o', '{}/sim_{}_'.format(self.getArg('outd'),seq)])
        self.bashout.write(artparam+'\n')
        os.system(artparam)   
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

if __name__ == "__main__":
  ttr=TreeToReads()

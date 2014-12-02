#!/usr/bin/python
import dendropy
from subprocess import call
import random
import os
import sys
import argparse
import numpy


class TreeToReads:
  def __init__(self, configfi='seqsim.cfg', run=1):
        """A method to read a tree, resolve polytomes, generate mutations and simulate reads."""
        self.configfi=configfi
        if run:
          self.runSims()
  def runSims(self):
            self.bashout = open('analysis.sh','w')
            self.readArgs()
            self.setDefaults()
            self.checkArgs()
            self.makeOut()
            self.readTree()
            self.readGenome()
            self.generateVarsites()
            self.readVarsites()
            self.selectMutsites()
#            self.mutGenomes()
            self.bashout.close()          
  def readArgs(self):
    if len(sys.argv) > 1:
      conf=sys.argv[1]
      try:
        config=open(conf).readlines()
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
  def makeOut(self):
      if not os.path.isdir(self.getArg('outd')):
          os.mkdir(self.getArg('outd'))
  def getArg(self,nam):
    return(self.config[self.argdict[nam]])
  def setDefaults(self): #TEMPORARY
            self.clustering=1
            self.clustPerc=0.25
            self.clustVar=50
            self.clustLoc=75
            self.alpha=2
            self.beta=2
  def checkArgs(self):
    self.argdict={'treepath':'treefile_path', 'nsnp':'number_of_snps', 'anchor_name':'anchor_name', 'genome':'anchor_genome_path', 'ratmat':'rate_matrix', 'freqmat':'freq_matrix', 'shape':'shape', 'errmod1':'error_model1', 'errmod2':'error_model2', 'cov':'coverage', 'outd':'output_dir'}
    for arg in self.argdict:
      if self.argdict[arg] not in self.config:
        print("{} is missing from the config file".format(self.argdict[arg]))
        sys.exit()
    try:
      int(self.getArg('nsnp'))
    except:
        print("number of SNPs {} could not be coerced to an integer".format(self.getArg('nsnp')))
        sys.exit()
    if not len(self.getArg('ratmat').split(','))==6:
        print("{} values in rate matrix, should be 6".format(len(self.getArg('ratmat').split(','))))
        sys.exit()
    if not len(self.getArg('freqmat').split(','))==4:
        print("{} values in freq matrix, should be 4".format(len(self.getArg('freqmat').split(','))))
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

  def readTree(self):
    #import tree from path
    taxa=dendropy.TaxonSet()
    tree=dendropy.Tree.get_from_path(self.getArg('treepath'),'newick',taxon_set=taxa)
    #fix polytomies
    tree.resolve_polytomies()
    self.outtree="{}/simtree.tre".format(self.getArg('outd'))
    tree.write(open(self.outtree,'w'),'newick',suppress_internal_node_labels=True)
    linrun="sed -i -e's/\[&U\]//' {}".format(self.outtree)
    self.bashout.write(linrun+'\n')
    os.system(linrun)

  def readGenome(self):
    genfas=open(self.getArg('genome')).readlines()
    crop=[lin[:-1] for lin in genfas[1:]]
    self.gen="".join(crop)
    self.genlen=len(self.gen)

  def generateVarsites(self):
## TODO make model variable
    self.simloc="{}/seqs_sim.txt".format(self.getArg('outd'))
    seqcall=" ".join(['seq-gen', '-l1000', '-n1', '-mGTR', '-a{}'.format(self.getArg('shape')), '-r{}'.format(self.getArg('ratmat')), '-f{}'.format(self.getArg('freqmat')), '<', '{}'.format(self.outtree), '>', '{}'.format(self.simloc)])
    os.system(seqcall)
    self.bashout.write(seqcall +'\n')

  def readVarsites(self):
    sims=open(self.simloc).readlines()[1:]
    bases=sims[0][10:]
    nucsets={}
    for i in range(len(bases)):
        nucsets[i]=set()
    for lin in sims:
      seq=lin[:10]
      bases=lin[10:-1]
      for i, nuc in enumerate(bases):
        nucsets[i].add(nuc)
    varsites=[]
    for i in nucsets:
      if len(nucsets[i]) == 2: #THIS LIMITS TO NO MULT HITS!! DO I want that?
        varsites.append(i)
    for vi in varsites:
        nucs=set()
    for lin in sims:
        bases=lin[10:]
        nuc=bases[vi]
        nucs.add(nuc)
    if len(nucs) == 0:
        print vi
    simseqs={} #What is happening here again?!
    for lin in sims:
      seq=lin[:10].strip()
      simseqs[seq]=[]
      bases=lin[10:-1]
      for i, nuc in  enumerate(bases):
        if i in varsites:
            simseqs[seq].append(nuc)
    self.seqnames=simseqs.keys()
    ref=simseqs[self.getArg('anchor_name')]
    print("Checkpoint 1")
    print(len(ref))
    self.sitepatts={}
    for nuc in ['A','G','T','C']:
	     self.sitepatts[nuc]=[]
    for i, nuc in enumerate(ref[:-1]):
      site={}
      nucs=set()
      for srr in simseqs:
          site[srr]=simseqs[srr][i]
          nucs.add(simseqs[srr][i])
      assert(len(nucs)>1)
      self.sitepatts[nuc].append(site)


  def selectMutsites(self):
    self.mutsite="{}/mutsites.txt".format(self.getArg('outd'))
    fi=open(self.mutsite,"w")
    rands=set()
    i=0
    while i < int(self.getArg('nsnp')): #This is the number of SNPs
      if i%100==0:
        print('.')
      if self.clustering and i < (self.genlen*self.clustPerc):
            i+=2
            ran=random.choice(range(self.genlen))
            rands.add(ran)
            fi.write(str(ran)+'\n')
            diff=int(random.gammavariate(self.alpha, self.beta)*20)
            if (random.choice([0, 1]) or (ran-diff < 0)) and (ran+diff < self.genlen):
              ranpair=ran+diff
            else:
              ranpair=ran-diff #Can run off of the end of the genome...
            rands.add(ranpair)
            fi.write(str(ranpair)+'\n')

      else:
            i+=1
            ran=random.choice(range(self.genlen))
            rands.add(ran)
            fi.write(str(ran)+'\n')
    self.mutlocs=rands
    # set of random numbers determining where the SNPs really fall. can be drawn from SNPlist or just random.
    fi.close()


# generate mutated genomes:
  def mutGenomes(self):
    self.mut_genos={}
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
        for nuc in self.gen:
                if ii%80==0:
                   genout.write('\n')
                ii+=1
                if ii in self.mutlocs:
                   patnuc[nuc]+=1
    #               print(nuc,patnuc[nuc])
                   patt=self.sitepatts[nuc][patnuc[nuc]]
                   genout.write(patt[seq])
                   self.mut_genos[seq].append(patt[seq])
    #               genout.write('X')
    #               print('MUT')
                else:
                    genout.write(nuc)
        genout.close()    

    def runART(self):
     #   print(' '.join(['art_illumina', '-1', '../simB/fullprofR1.txt', '-2', '../simB/fullprofR2.txt', '-p', '-sam', '-i', 'sim_{}.fasta'.format(seq), '-l', '150', '-f', '20', '-m', '350', '-s', '130', '-o', 'sim_{}_'.format(seq)]))
     #   call(['art_illumina', '-1', args['error_model1'], '-2', args['error_model2'], '-p', '-sam', '-i', 'sim_{}.fasta'.format(seq), '-l', '150', '-f', args['coverage'], '-m', '350', '-s', '130', '-o', 'sim_{}_'.format(seq)])
        artparam=' '.join(['art_illumina', '-1', self.config[self.argdict[errmod2]], '-2', self.config[self.argdict[errmod2]], '-p', '-sam', '-i', 'sim_{}.fasta'.format(seq), '-l', '150', '-f', self.config[self.argdict[cov]], '-m', '350', '-s', '130', '-o', '{}/sim_{}_'.format(self.getArg('outd'),seq)])
        os.system(artparam)   
    '''
    #'-l', '150', <- average read length
    # '-f', '20', <-need to get fold coverage
    #'-m', '350', <- fragment size?? need to check
    #'-s', '130', <- stdev fragment size...
     
    # pull a random number fromt the length of the genome,
    #art_illumina -1 fullprof.txt -2 fullprof.txt -p -sam -i SIM_.fasta -l 150 -f 20 -m 2050 -s 50 -o sim_test
    #Map those mutations back onto the genomes....
    #Get distributaion of SNP locations from real data?!!
    # simulate reads for each of those genomes.
    '''


ttr=TreeToReads()

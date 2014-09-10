import dendropy
from subprocess import call
import random
import os
import sys

import argparse
if len(sys.argv) > 1:
    conf=sys.argv[1]
    try:
      config=open(conf).readlines()
    except:
      print("couldn't open config file")
      sys.exit()
else:
  try:
    config=open("seqsim.cfg").readlines()
  except:
      print("config file seqsim.cfg not found")
      sys.exit()


args={}
for lin in config:
  lii=lin.split('=')
  args[lii[0].strip()]=lii[-1].split('#')[0].strip()

treepath='treefile_path'
nsnp='number_of_snps'
anchor_name='anchor_name'
genome='anchor_genome_path'
ratmat='rate_matrix'
freqmat='freq_matrix'
shape='shape'

argslist=['treefile_path', 'number_of_snps', 'anchor_name', 'anchor_genome_path', 'rate_matrix', 'freq_matrix', 'shape', 'error_model1']
for item in argslist:
  if item not in args:
    print("{} is missing from the config file".format(item))
    sys.exit()

try:
  int(args[nsnp])
except:
  print("number of SNPs {} could not be coerced to an integer".format(args[nsnp]))
  sys.exit()

if not len(args[ratmat].split(','))==6:
  print("{} values in rate matrix, should be 6".format(len(args[ratmat].split(','))))
  sys.exit()

if not len(args[freqmat].split(','))==4:
  print("{} values in freq matrix, should be 4".format(len(args[freqmat].split(','))))
  sys.exit()  

try:
  float(args[shape])
except:
  print("shape parameter {} could not be coerced to a float".format(args[shape]))
  sys.exit() 

try:
  open(args[treepath])
except:
  print("could not open treefile {}".format(args[treepath]))
  sys.exit() 

try:
  open(args[genome])
except:
  print("could not open anchor genome {}".format(args[genome]))
  sys.exit() 


#import tree from path
taxa=dendropy.TaxonSet()
tree=dendropy.Tree.get_from_path(args[treepath],'newick',taxon_set=taxa)

#fix polytomies
tree.resolve_polytomies()
tree.write(open("simtree.tre",'w'),'newick',suppress_internal_node_labels=True)
os.system("sed -i -e's/\[&U\]//' simtree.tre")


genfas=open(args[genome]).readlines()
crop=[lin[:-1] for lin in genfas[1:]]
gen="".join(crop)
genlen=len(gen)
print args[shape]

## TODO make model variable
seqcall=" ".join(['seq-gen', '-l1000', '-n1', '-mGTR', '-a{}'.format(args[shape]), '-r{}'.format(args[ratmat]), '-f{}'.format(args[freqmat]), '<', 'simtree.tre', '>', 'seqs_sim.txt'])
os.system(seqcall)

sims=open("seqs_sim.txt").readlines()[1:]

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
  if len(nucsets[i]) == 2:
    varsites.append(i)

for vi in varsites:
  nucs=set()
  for lin in sims:
    bases=lin[10:]
    nuc=bases[vi]
    nucs.add(nuc)
  if len(nucs) == 0:
    print vi


genos={}
for lin in sims:
    seq=lin[:10].strip()
    genos[seq]=[]
    bases=lin[10:-1]
    for i, nuc in  enumerate(bases):
        if i in varsites:
            genos[seq].append(nuc)

 
ref=genos[args[anchor_name]]

sitepatts={}
for nuc in ['A','G','T','C']:
	sitepatts[nuc]=[]


for i, nuc in enumerate(ref[:-1]):
    site={}
    nucs=set()
    for srr in genos:
        site[srr]=genos[srr][i]
        nucs.add(genos[srr][i])
    assert(len(nucs)>1)
    sitepatts[nuc].append(site)


fi=open("mutlocs.txt","w")
rands=set()
i=0


while i < int(args[nsnp]): #This is the number of SNPs
    i+=1
    ran=random.choice(range(genlen))
    rands.add(ran)
    fi.write(str(ran)+'\n')
# set of random numbers determining where the SNPs really fall. can be drawn from SNPlist or just random.
fi.close()


# generate mutated genomes:

for seq in genos:
    genout=open("sim_{}.fasta".format(seq),'w')
    patnuc={}
    patnuc['A']=0
    patnuc['G']=0
    patnuc['T']=0
    patnuc['C']=0
    ii=0
    genout.write(">SIM_{}\n".format(seq))
    for lin in genome[1:]:
        for nuc in lin[:-1]:
            ii+=1
            if ii in rands:
               patnuc[nuc]+=1
#               print(nuc,patnuc[nuc])
               patt=sitepatts[nuc][patnuc[nuc]]
               genout.write(patt[seq])
#               genout.write('X')
#               print('MUT')
            else:
                genout.write(nuc)
        genout.write('\n')
    genout.close()    
 #   print(' '.join(['art_illumina', '-1', '../simB/fullprofR1.txt', '-2', '../simB/fullprofR2.txt', '-p', '-sam', '-i', 'sim_{}.fasta'.format(seq), '-l', '150', '-f', '20', '-m', '350', '-s', '130', '-o', 'sim_{}_'.format(seq)]))
    call(['art_illumina', '-1', '../simA/fullprofR1.txt', '-2', '../simA/fullprofR2.txt', '-p', '-sam', '-i', 'sim_{}.fasta'.format(seq), '-l', '150', '-f', '45', '-m', '350', '-s', '130', '-o', 'sim_{}_'.format(seq)])
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

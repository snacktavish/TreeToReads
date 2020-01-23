#!/usr/bin/env python

import dendropy

import sys
ref = sys.argv[1]
sim = sys.argv[2]

inputseq = dendropy.DnaCharacterMatrix.get(file=open(ref), schema="fasta")
outputseq = dendropy.DnaCharacterMatrix.get(file=open(sim), schema="fasta")

assert len(inputseq) == len(outputseq)
i = 0
for item in inputseq: 
  assert(inputseq[i].symbols_as_string() == outputseq[i].symbols_as_string().replace('-',''))
  i+=1

print("check refgenome test passed\n")
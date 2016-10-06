#!/bin/bash

echo "Running basic tests"

python treetoreads.py tests/multi.config >> tests/test.log

python  tests/check_refgenome.py tests/input/mini_ref_multi.fasta tests/multi/fasta_files/sim_gi.fasta 

rm -r tests/multi

python treetoreads.py tests/multi_indel.config >> tests/test.log

python  tests/check_refgenome.py tests/input/mini_ref_multi.fasta tests/multi_indel/fasta_files/sim_gi.fasta

rm -r tests/multi_indel

echo "integration"
python tests/integration_test.py 

echo "tree reconstruction from fasta"
python tests/tree_from_fasta.py


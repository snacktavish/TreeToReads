echo "runnuig basic tests"

echo "tree reconstrutcion from fasta"
python tests/tree_from_fasta.py
echo "simple clustering run"
python tests/clustering.py
echo "tree reconstrutcion via SNP pipeline"
python tests/tree_reconstruction_test.py

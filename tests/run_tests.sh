echo "runnuig basic tests"

echo "tree reconstrutcion from fasta"
python tree_from_fasta.py
echo "simple clustering run"
python clustering.py
echo "tree reconstrutcion via SNP pipeline"
python tree_reconstruction_test.py
echo "Running basic tests"

echo "tree reconstruction from fasta"
python tests/tree_from_fasta.py
echo "simple clustering run"
python tests/clustering.py
echo "integration"
python tests/integration_test.py

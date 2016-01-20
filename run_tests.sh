echo "Running basic tests"

echo "tree reconstruction from fasta"
python tests/tree_from_fasta.py
echo "integration"
python tests/integration_test.py

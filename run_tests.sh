echo "Running basic tests"

echo "integration"
python tests/integration_test.py

echo "tree reconstruction from fasta"
python tests/tree_from_fasta.py


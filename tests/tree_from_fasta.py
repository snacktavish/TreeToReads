import unittest
from treetoreads import TreeToReads
import subprocess
import dendropy
import subprocess
import os


def compare_trees(expected,estimated):
 #   assert(estimated.euclidean_distance(expected)<= branch lengths are hard to test. TODO figure out how.
    taxon_namespace = dendropy.TaxonSet()
    exp_tree = dendropy.Tree.get_from_path(
        expected,
        "newick",
        taxon_set=taxon_namespace)
    est_tree = dendropy.Tree.get_from_path(
        estimated,
        "nexus",
        taxon_set=taxon_namespace)
    return(est_tree.symmetric_difference(exp_tree))

if __name__ == '__main__':
    ttr=TreeToReads(configfi='/home/ejmctavish/Documents/FDA/TreetoReads/tests/input/tree_from_fasta.cfg',run=0)
    print("outd is {}".format(ttr.outd))
    ttr.outd="/home/ejmctavish/Documents/FDA/TreetoReads/tests/tree_from_fasta"
    ttr.makeOut()
    ttr.mutGenomes()
    print(ttr.getArg('genome'))
    # enter the directory like this:
    os.system("cat {}/*.fasta > {}/mini.aln".format(ttr.outd, ttr.outd))
    os.system("Garli tests/garli_mini.conf")
    assert(compare_trees('tests/input/simple_expected.tre','tests/tree_from_fasta/tff.best.tre')==0)
    assert(compare_trees('tests/input/alt.tre','tests/tree_from_fasta/tff.best.tre')!=0)
    print('Topology is correct')
    os.system("rm -r /home/ejmctavish/Documents/FDA/TreetoReads/tests/tree_from_fasta")

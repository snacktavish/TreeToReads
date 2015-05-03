from treetoreads import TreeToReads
import subprocess

def runSNPpipe(ttr):
    print(" ".join(['run_snppipe.sh',ttr.outd, ttr.getArg('genome')]))
    subprocess.Popen(['./run_snppipe.sh',ttr.outd, ttr.getArg('genome')])

if __name__ == '__main__':
    ttr=TreeToReads(configfi='tests/input/simple.cfg',run=0)
    #ttr._vargen=1
    #ttr.simloc="tests/input/test_seqs_sim.txt"
    ttr.readTree()
#    print(ttr.seqnames)
    ttr.readVarsites()
 #   print(ttr.sitepatts)
    ttr.mutGenomes()
    ttr.runART()
 #   print("outd is {}".format(ttr.outd))
 #   print(ttr.getArg('genome'))
 #   runSNPpipe(ttr)
     cat sim_CFSAN000189.fasta <(echo -e '\n') sim_CFSAN000211.fasta <(echo -e '\n') sim_CFSAN000228.fasta  <(echo -e '\n') sim_CFSAN000958.fasta  <(echo -e '\n') sim_CFSAN000961.fasta <(echo -e '\n') sim_CFSAN000191.fasta <(echo -e '\n') sim_CFSAN000212.fasta <(echo -e '\n') sim_CFSAN000669.fasta <(echo -e '\n') sim_CFSAN000960.fasta <(echo -e '\n') sim_CFSAN001140.fasta > test.aln
cat sim_A.fasta <(echo -e '\n') sim_B.fasta <(echo -e '\n') sim_C.fasta  <(echo -e '\n') sim_D.fasta  <(echo -e '\n') sim_ref.fasta < mini.aln
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
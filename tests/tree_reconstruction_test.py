import unittest
from treetoreads import TreeToReads
import subprocess

def runSNPpipe(ttr):
    print(" ".join(['run_snppipe.sh',ttr.outd, ttr.getArg('genome')]))
    subprocess.Popen(['./run_snppipe.sh',ttr.outd, ttr.getArg('genome')])


if __name__ == '__main__':
    ttr=TreeToReads(configfi='tests/input/simple.cfg',run=0)
    ttr.runART()
    print("outd is {}".format(ttr.outd))
    print(ttr.getArg('genome'))
    runSNPpipe(ttr)

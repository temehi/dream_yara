#!/usr/bin/env python2

import sys
import os
import shutil

dirName = 'input/'+ sys.argv[2] + '-binned-genomes'
if sys.argv[1] == 'split':
    # ============================================================
    # Split the genomes in to 64 bins.
    # ============================================================
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    binFastaLen = 101 # each fasta entry is 100 lines long (101 including the header)
    inputFasta = open('input/'+ sys.argv[2] + '-genomes.fa', 'r').read().split('\n')
    for b in range(64):
        # First, get the slice
        binFasta = inputFasta[b*binFastaLen:(b+1)*binFastaLen]
        # Now open one file per bin and dump the part
        outputFasta = open(dirName + '/' + str(b) + '.fa', 'w')
        outputFasta.write('\n'.join(binFasta))
        outputFasta.close()
elif sys.argv[1] == 'clear':
    # ============================================================
    # remove spllited files
    # ============================================================
    if os.path.exists(dirName):
        shutil.rmtree(dirName)
else:
    print ("Unknown operation " + sys.argv[1])
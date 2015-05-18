#!/usr/bin/env python

from argparse import ArgumentParser, FileType

import numpy as np
from numpy.random import exponential, choice, uniform
from scipy.special import binom

class Lineage:
    def __init__(self, sequence, parent=None, height=float("inf")):
        self.sequence = sequence
        self.parent = parent
        self.height = height

        self.subs = []
        self.maxHeight = 0
        self.ancSeq = sequence[:]

    def sim(self, height):
        """Impute ancestral sequences."""
        
        if self.parent == None:
            t = self.maxHeight
            while True:
                t += exponential(1.0/len(self.sequence))
                if (t > height):
                    break

                i = choice(len(self.sequence))
                self.ancSeq[i] = (self.ancSeq[i] + choice(range(1,4))) % 4
                self.subs.append((t, i, self.ancSeq[i]))

        self.maxHeight = height


    def getHD(self, height, otherSequence):

        if (height > self.maxHeight):
            self.sim(height)
            thisSeq = self.ancSeq
        else:
            thisSeq = self.sequence[:]
            for subHeight, subIdx, sub in self.subs:
                if subHeight > height:
                    break

                thisSeq[subIdx] = sub

        return sum(map(lambda x: x[0]!=x[1], zip(otherSequence, thisSeq)))

def probHD(h, t, L):
    """Probability of evolving to sequence h away from start in time t."""
    p = 0.75*(1.0-np.e**(-t*4.0/3.0))

    #return binom(L,h)*((1-p)**(L-h))*(p**h)
    return ((1-p)**(L-h))*(p**h)


def readFASTA(inFile):
    seqs = {}
    alphabet = 'GCTA'

    header = None
    for line in inFile:
        if line.startswith('>'):
            header = line[1:].strip()
        else:
            seqs[header] = map(lambda c: alphabet.index(c), list(line.strip()))

    return seqs


## MAIN ##

if __name__ == '__main__':

    parser = ArgumentParser(description='Direct tree sampler.')
    parser.add_argument("fasta", type=FileType('r'), help='FASTA file (contemporaneously sampled)')

    args = parser.parse_args()

    seqs = readFASTA(args.fasta)

    lineage0 = Lineage(seqs.values()[0])
    print lineage0.getHD(0, seqs.values()[1])

    L = len(seqs.values()[0])

    T = None
    while True:
        T = exponential(0.01)
        h = lineage0.getHD(T, seqs.values()[1])
        p = probHD(h, T, L)/1e-70
        print p
        if uniform() < p:
            break

    print h,T

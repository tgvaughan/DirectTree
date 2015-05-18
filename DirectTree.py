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

    L = len(seqs.values()[0])

    print "T"
    for s in range(1000):
        lineage0 = Lineage(seqs.values()[0])
        T = 0.1
        h = lineage0.getHD(T, seqs.values()[1])
        P = probHD(h, T, L)*np.e**(-T)

        for iter in range(100):
            #Tp = T + uniform(-0.001,0.001)
            f = uniform(0.8, 1.0/0.8)
            Tp = T*f

            if Tp<0:
                alpha = 0
            else:
                hp = lineage0.getHD(Tp, seqs.values()[1])
                Pp = probHD(hp, Tp, L)*np.e**(-Tp)
                alpha = min(Pp/P/f, 1.0)

            if alpha == 1 or (alpha > 0 and uniform()<alpha):
                T = Tp
                P = Pp

        print T

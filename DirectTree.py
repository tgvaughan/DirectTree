#!/usr/bin/env python

from argparse import ArgumentParser, FileType
from sys import stdout

import numpy as np
from numpy.random import exponential, uniform, randint
from scipy.special import binom

from matplotlib import pyplot as plt

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

                i = randint(len(self.sequence))
                self.ancSeq[i] = (self.ancSeq[i] + randint(1,4)) % 4
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

def logProbHD(h, t, L):
    """Probability of evolving to sequence h away from start in time t."""
    p = 0.25*(1.0-np.e**(-t*4.0/3.0))

    return (L-h)*np.log(1-3*p) + h*np.log(p)


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
    parser.add_argument("--theta", type=float, default=1.0, help='Parameter value.')
    parser.add_argument("-o", dest="outfile", type=FileType('w'), default=stdout, help='Output file (default stdout)')

    args = parser.parse_args()

    seqs = readFASTA(args.fasta)

    L = len(seqs.values()[0])

    args.outfile.write("T\n")
    for s in range(1000):
        T = args.theta
        lineage0 = Lineage(seqs.values()[0])
        h = lineage0.getHD(T, seqs.values()[1])
        logP = logProbHD(h, T, L) - T/args.theta

        #Tvec = [T]

        for iter in range(200):

            if (uniform() < 0.5):
                f = uniform(0.5, 1.0/0.5)
            else:
                f = uniform(0.9, 1.0/0.9)

            Tp = T*f

            if Tp<0:
                logalpha = float('-inf')
            else:
                hp = lineage0.getHD(Tp, seqs.values()[1])
                logPp = logProbHD(hp, Tp, L) - Tp/args.theta
                logalpha = logPp - logP - np.log(f)

            if logalpha>0  or uniform()<np.e**logalpha:
                T = Tp
                logP = logPp
    
            #Tvec.append(T)
        args.outfile.write(str(T) + "\n")
        args.outfile.flush()

    #plt.plot(range(len(Tvec)), Tvec)
    #plt.show()

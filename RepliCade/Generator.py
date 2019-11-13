'''

Python File that contains a Generator class for
generating DNA sequences with various different generator
functions and parameters

'''

import random

class Generator:

    def __init__(self):
        self.sequence = ""


    def generateSequence(self,seqLen):
        '''Function that generates a random DNA sequence of a user specified length'''
        nucleotides = ['A','G','T','C']
        sequence = ""
        for i in range(0, seqLen):
            sequence+=random.choice(nucleotides)
        self.sequence = sequence
        return sequence


    def generateSequenceGCContent(self, seqLen, GCContent):
        '''Function that generates a random DNA sequence of a user specified length
           and user specified GC content'''
        if (GCContent > 1 or GCContent < 0):
            return "Invalid GC Content. GC Content cannot be more than 100 percent or less than 0 percent"
        gcAr = ['G','C']
        other = ['A','T']
        ar = []
        gc = seqLen*GCContent
        for i in range(0, seqLen):
            if i < gc:
                ar.append(random.choice(gcAr))
            else:
                ar.append(random.choice(other))
        random.shuffle(ar)
        self.sequence = ''.join(ar)
        return self.sequence


    def getSequence(self):
        return self.sequence






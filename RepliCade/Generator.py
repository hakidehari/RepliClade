'''

Python File that contains a Generator class for
generating DNA sequences with various different generator
functions and parameters

'''

import random

class Generator:

    def __init__(self):
        self.sequence = ""
        self.nucleotides = ['A','G','T','C']

    '''Function that generates a random DNA sequence of a user specified length'''
    def generateSequence(self,seqLen):
        sequence = ""
        for i in range(0, seqLen):
            sequence+=self.nucleotides[random.randint(0, 3)]

        self.sequence = sequence
        return sequence


    '''Function that generates a random DNA sequence of a user specified length
        and user specified GC content'''
    def generateSequenceGCContent(self, seqLen, GCContent):
        if (GCContent > 1 or GCContent < 0):
            return "Invalid GC Content. GC Content cannot be more than 100 percent or less than 0 percent"
        ar = []
        gc = seqLen*GCContent
        for i in range(0, seqLen):
            r = random.randint(0, 1)
            if i < gc:
                if r == 0:
                    ar.append('G')
                else:
                    ar.append('C')
            else:
                if r == 0:
                    ar.append('A')
                else:
                    ar.append('T')
        random.shuffle(ar)
        self.sequence = ''.join(ar)
        return self.sequence

    def getSequence(self):
        return self.sequence






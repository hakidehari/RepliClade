
import random


class Sequence:

    def __init__(self, sequence):
        self.sequence = sequence



    '''Modifies a sequence with a random nucleotide substitution at the user specified
        position'''
    def modify(self, position):
        if position > len(self.sequence) - 1 or position < 0:
            return "Index provided is out of bounds"
        nucleotides = ['A','G','C','T']
        ar = list(self.sequence)
        ar[position] = nucleotides[random.randint(0, 3)]
        self.sequence = ''.join(ar)
        #return self.sequence

    '''Modifies a sequence with a random nucleotide substitution at a random position'''
    def randomModify(self):
        nucleotides = ['A', 'G', 'C', 'T']
        ar = list(self.sequence)
        ar[random.randint(0, len(self.sequence)-1)] = nucleotides[random.randint(0, 3)]
        self.sequence = ''.join(ar)
        return self.sequence

    def getSequence(self):
        return self.sequence




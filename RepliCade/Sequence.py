
import random
from Generator import Generator

"""

This class provides the structure for DNA sequences
The class contains certain functions that modify the sequence if required
More detailed functions for indels and other modifications will be added in the future

"""

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
        return self.sequence


    '''Modifies a sequence with a random nucleotide substitution at a random position'''
    def randomModify(self):
        nucleotides = ['A', 'G', 'C', 'T']
        ar = list(self.sequence)
        ar[random.randint(0, len(self.sequence)-1)] = random.choice(nucleotides)
        self.sequence = ''.join(ar)
        return self.sequence


    '''
        Generate indels to insert of random lengths (for right now, I've chosen
        lengths between 2 and 20) using probabilities of length insertion
        More complex probability models will be added at a later time
    '''
    def insertIndel(self):
        indelLength = random.randint(2, 20)
        indel = Generator.generateSequence(Generator, indelLength)
        indelPosition = random.randint(0, len(self.sequence))
        self.sequence = self.sequence[:indelPosition] + indel + self.sequence[indelPosition:]
        return self.sequence

    def getSequence(self):
        return self.sequence




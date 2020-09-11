"""

This class provides the structure for DNA sequences
The class contains certain functions that modify the sequence if required
More detailed functions for indels and other modifications will be added in the future

"""

import random
from Generator import Generator

class Sequence:

    def __init__(self, sequence):
        self.sequence = sequence
        self.generator = Generator()


    def modify(self, position):
        '''Modifies a sequence with a random nucleotide substitution at the user specified
        position'''
        if position > len(self.sequence) - 1 or position < 0:
            return "Index provided is out of bounds"
        nucleotides = ['A','G','C','T']
        ar = list(self.sequence)
        ar[position] = nucleotides[random.randint(0, 3)]
        self.sequence = ''.join(ar)


    def random_modify(self):
        '''Modifies a sequence with a random nucleotide substitution at a random position'''
        nucleotides = ['A', 'G', 'C', 'T']
        ar = list(self.sequence)
        ar[random.randint(0, len(self.sequence)-1)] = random.choice(nucleotides)
        self.sequence = ''.join(ar)


    
    def insert_indel(self):
        '''
        Generate indels to insert of random lengths (for right now, I've chosen
        lengths between 2 and 20) using probabilities of length insertion
        More complex probability models will be added at a later time
        
        '''
        indel_length = random.randint(2, 20)
        indel = self.generator.generate_sequence(indel_length)
        indel_position = random.randint(0, len(self.sequence))
        self.sequence = self.sequence[:indel_position] + indel + self.sequence[indel_position:]
        return self.sequence


    def get_sequence(self):
        '''Returns sequence of sequence object'''
        return self.sequence



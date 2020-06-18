'''

Python File that contains a Generator class for
generating DNA sequences with various different generator
functions and parameters

'''

import random

class Generator:

    def __init__(self):
        self.sequence = ""


    def generate_sequence(self,seqLen):
        '''Function that generates a random DNA sequence of a user specified length'''
        nucleotides = ['A','G','T','C']
        sequence = ""
        for i in range(0, seqLen):
            sequence+=random.choice(nucleotides)
        self.sequence = sequence
        return sequence


    def generate_sequence_gc_content(self, seq_len, gc_content):
        '''Function that generates a random DNA sequence of a user specified length
           and user specified GC content'''
        if (gc_content > 1 or gc_content < 0):
            return "Invalid GC Content. GC Content cannot be more than 100 percent or less than 0 percent"
        gc_arr = ['G','C']
        other = ['A','T']
        ar = []
        gc = seq_len*gc_content
        for i in range(0, seqLen):
            if i < gc:
                ar.append(random.choice(gc_arr))
            else:
                ar.append(random.choice(other))
        random.shuffle(ar)
        self.sequence = ''.join(ar)
        return self.sequence


    def get_sequence(self):
        return self.sequence






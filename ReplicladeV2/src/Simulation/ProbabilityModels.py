"""
A modular file with Probability Model functions for many different 
genome types/scenarios
"""
import random

nucleotides = ['A','G','C','T']

def mutate_sequence_random(sequence):
    '''
        A random mutation using random probabilities of a nucleotide sequence
    '''
    s = list(sequence)
    for i in range(0, len(s)):
        threshold = random.random()
        random_val = random.random()
        if random_val < threshold:
            s[i] = random.choice(nucleotides)
    return ''.join(s)
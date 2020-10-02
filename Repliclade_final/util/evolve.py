import random
import numpy as np
import math

class JukesCantor(object):


    def __init__(self):
        '''
        2D probability matrix for nucleotide substitutions
        in the Jukes and Cantor model, all of the probabilities from one nucleotide to the other
        are .25
        '''
        self.alpha = .001
        self.t = 0

        self.calculate_matrix(self.alpha, self.t)

        self.seq_list = ['A', 'T', 'C', 'G']


    def calculate_matrix(self, alpha, t):
        '''
        Markov model which defines the probability substitution matrix
        in the current generation
        '''
        same_nuc = .25 + .75*(math.e**(-4*alpha*t))
        diff_nuc = .25-.25*(math.e**(-4*alpha*t))
        self.prb_matrix =  {
            #     A    T    C    G
            'A': [same_nuc, diff_nuc, diff_nuc, diff_nuc],
            'T': [diff_nuc, same_nuc, diff_nuc, diff_nuc],
            'C': [diff_nuc, diff_nuc, same_nuc, diff_nuc],
            'G': [diff_nuc, diff_nuc, diff_nuc, same_nuc]
        }


    def evolve(self, seq):
        '''
        Takes an input sequence and uses the Jukes and Cantor model
        to evolve the sequence
        '''

        ret_seq = ''
        for i in range(len(seq)):
            cur = seq[i]
            if cur == '-':
                ret_seq += cur
                continue
            for j in range(len(self.prb_matrix[cur])):
                if self.prb_matrix[cur][j] != 0:
                    roll = random.random()
                    if roll <= self.prb_matrix[cur][j]:
                        cur = self.seq_list[j]
                        break
            ret_seq += cur
        self.t += 1
        self.calculate_matrix(self.alpha, self.t)
        return ret_seq


class Kimura(object):

    
    def __init__(self):
        '''
        2D probability matrix for nucleotide substitutions
        in the Jukes and Cantor model, all of the probabilities from one nucleotide to the other
        are .25
        '''
        self.alpha = .001
        self.t = 0

        self.calculate_matrix(self.alpha, self.t)

        self.seq_list = ['A', 'T', 'C', 'G']


    def calculate_matrix(self, alpha, t):
        '''
        Markov model which defines the probability substitution matrix
        in the current generation
        '''
        same_nuc = .25 + .75*(math.e**(-4*alpha*t))
        diff_nuc = .25-.25*(math.e**(-4*alpha*t))
        self.prb_matrix =  {
            #     A    T    C    G
            'A': [same_nuc, diff_nuc/2, diff_nuc/2, diff_nuc],
            'T': [diff_nuc/2, same_nuc, diff_nuc, diff_nuc/2],
            'C': [diff_nuc/2, diff_nuc, same_nuc, diff_nuc/2],
            'G': [diff_nuc, diff_nuc/2, diff_nuc/2, same_nuc]
        }
        
    
    def evolve(self, seq):
        '''
        Takes an input sequence and uses the Kimura model
        to evolve the sequence
        '''

        ret_seq = ''
        for i in range(len(seq)):
            cur = seq[i]
            if cur == '-':
                ret_seq += cur
                continue
            for j in range(len(self.prb_matrix[cur])):
                if self.prb_matrix[cur][j] != 0:
                    roll = random.random()
                    if roll <= self.prb_matrix[cur][j]:
                        cur = self.seq_list[j]
                        break
            ret_seq += cur
        self.t += 1
        self.calculate_matrix(self.alpha, self.t)
        return ret_seq


class Blaisdell(object):


    def __init__(self):
        pass
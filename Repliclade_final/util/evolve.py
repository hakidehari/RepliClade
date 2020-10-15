import random
import numpy as np
import math
from util.seq_util import SequenceUtil

seq_util = SequenceUtil()

class JukesCantor(object):


    def __init__(self):
        '''
        2D probability matrix for nucleotide substitutions
        in the Jukes and Cantor model, all of the probabilities from one nucleotide to the other
        are .25
        '''
        self.alpha = .0001
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

    
    def evolve_cr(self, seq, c_regions):
        ret_seq = ''
        for i in range(len(seq)):
            cur = seq[i]
            if seq_util.check_for_cr(i, c_regions):
                ret_seq += cur
                continue
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
        self.alpha = .0001
        self.t = 0
        self.beta= self.alpha / 3

        self.calculate_matrix(self.alpha, self.beta, self.t)

        self.seq_list = ['A', 'T', 'C', 'G']


    def calculate_matrix(self, alpha, beta, t):
        '''
        Markov model which defines the probability substitution matrix
        in the current generation
        '''
        transition = .25 + .25*(math.e**(-4*beta*t)) - .5*(math.e**(-2*(alpha + beta)*t))
        transversion = .25 - .25*(math.e**(-4*beta*t))
        same_nuc = 1 - transition - 2*transversion
        self.prb_matrix =  {
            #     A    T    C    G
            'A': [same_nuc, transversion, transversion, transition],
            'T': [transversion, same_nuc, transition, transversion],
            'C': [transversion, transition, same_nuc, transversion],
            'G': [transition, transversion, transversion, same_nuc]
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
        self.calculate_matrix(self.alpha, self.beta, self.t)
        return ret_seq

    
    def evolve_cr(self, seq, c_regions):
        '''
        Evolves the sequence while considering conserved regions
        '''
        ret_seq = ''
        for i in range(len(seq)):
            cur = seq[i]
            if seq_util.check_for_cr(i, c_regions):
                ret_seq += cur
                continue
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
        self.calculate_matrix(self.alpha, self.beta, self.t)
        return ret_seq


class Blaisdell(object):


    def __init__(self):
        pass


class Felsenstein(object):


    def __init__(self, seq):
        frequencies = seq_util.get_nuc_count(seq)
        self.A = float(frequencies['A'] / len(seq))
        self.C = float(frequencies['C'] / len(seq))
        self.G = float(frequencies['G'] / len(seq))
        self.T = float(frequencies['T'] / len(seq))
        self.t = 0
        self.calculate_matrix(self.t)
        self.seq_list = ['A', 'T', 'C', 'G']
        
        
    def calculate_matrix(self, t):
        '''
        Markov model which defines the probability substitution matrix
        in the current generation
        '''
        u = float(1 / (1 - self.A**2 - self.C**2 - self.G**2 - self.T**2))
        
        self.prb_matrix =  {
            #     A    T    C    G
            'A': [math.e**(-1*u*self.t) + self.A*(1 - math.e**(-1*u*self.t)), self.T*(1-math.e**(-1*u*self.t)), self.C*(1-math.e**(-1*u*self.t)), self.G*(1-math.e**(-1*u*self.t))],
            'T': [self.A*(1-math.e**(-1*u*self.t)), math.e**(-1*u*self.t) + self.T*(1 - math.e**(-1*u*self.t)), self.C*(1-math.e**(-1*u*self.t)), self.G*(1-math.e**(-1*u*self.t))],
            'C': [self.A*(1-math.e**(-1*u*self.t)), self.T*(1-math.e**(-1*u*self.t)), math.e**(-1*u*self.t) + self.C*(1 - math.e**(-1*u*self.t)), self.G*(1-math.e**(-1*u*self.t))],
            'G': [self.A*(1-math.e**(-1*u*self.t)), self.T*(1-math.e**(-1*u*self.t)), self.C*(1-math.e**(-1*u*self.t)), math.e**(-1*u*self.t) + self.G*(1 - math.e**(-1*u*self.t))]
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
        self.calculate_matrix(self.t)
        return ret_seq


    def evolve_cr(self, seq, c_regions):
        '''
        Evolves the sequence while considering conserved regions
        '''
        ret_seq = ''
        for i in range(len(seq)):
            cur = seq[i]
            if seq_util.check_for_cr(i, c_regions):
                ret_seq += cur
                continue
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
        self.calculate_matrix(self.t)
        return ret_seq
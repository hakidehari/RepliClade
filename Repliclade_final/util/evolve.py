import random

class JukesCantor(object):


    def __init__(self):
        #2D probability matrix for nucleotide substitutions
        #in the Jukes and Cantor model, all of the probabilities from one nucleotide to the other
        #are .25
        self.prb_matrix = {
            #     A    T    C    G
            'A': [0, .25, .25, .25],
            'T': [.25, 0, .25, .25],
            'C': [.25, .25, 0, .25],
            'G': [.25, .25, .25, 0]
        }

        self.seq_list = ['A', 'T', 'C', 'G']


    def evolve(self, seq):
        '''
        Takes an input sequence and uses the Jukes and Cantor model
        to evolve the sequence
        '''

        ret_seq = ''
        for i in range(len(seq)):
            cur = seq[i]
            for j in range(len(self.prb_matrix[cur])):
                if self.prb_matrix[cur][j] != 0:
                    roll = random.random()
                    if roll <= self.prb_matrix[cur][j]:
                        cur = self.seq_list[j]
            ret_seq += cur
        return ret_seq


class Kimura(object):

    
    def __init__(self):
        self.prb_matrix = {
            #     A    T    C    G
            'A': [0, .125, .125, .25],
            'T': [.125, 0, .25, .125],
            'C': [.125, .25, 0, .125],
            'G': [.25, .125, .125, 0]
        }

        self.seq_list = ['A', 'T', 'C', 'G']
        
    
    def evolve(self, seq):
        '''
        Takes an input sequence and uses the Kimura model
        to evolve the sequence
        '''

        ret_seq = ''
        for i in range(len(seq)):
            cur = seq[i]
            for j in range(len(self.prb_matrix[cur])):
                if self.prb_matrix[cur][j] != 0:
                    roll = random.random()
                    if roll <= self.prb_matrix[cur][j]:
                        cur = self.seq_list[j]
            ret_seq += cur
        return ret_seq


class Blaisdell(object):


    def __init__(self):
        pass
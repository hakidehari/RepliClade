"""

This is the Main file for Repliclade where the program will be run
All dependencies and classes will be compiled here and ran accordingly

"""

import random
from Generator import Generator
from Sequence import Sequence
from Connector import Connector
from EyeOh import EyeOh
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import Phylo
from ProbabilityModels import ProbabilityModels
import os


class Simulation:


    def performMSA(self):
        '''
                Performs multiple sequence alignment after simulation is run
                using clustalW2, one of the fastest MSA algorithms out there.
                Once finished, displays the sequences and their alignment
        '''
        print("Performing Multiple Sequence Alignment.......")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        clustalw_exe = os.path.join(dir_path, "clustalw2")
        clustalw_cline = ClustalwCommandline(clustalw_exe, infile="influenza.fasta")
        stdout, stderr = clustalw_cline()
        align = AlignIO.read("influenza.aln", "clustal")
        print(align)
        print("Alignment Complete.")



    def displayTree(self):
        '''Creates and displays a phylogenetic tree based on the Multiple Sequence Alignment scores'''
        tree = Phylo.read("influenza.dnd", "newick")
        tree.rooted = True
        Phylo.draw(tree)


    def getSequences(self, gene):
        '''Gathers sequences from GenBank'''
        con = Connector()
        seq = con.getGeneData(gene)
        return seq


    def printSequences(self, sequences):
        '''Prints the sequences returned from the simulation'''
        ar = []
        for sequence in sequences:
            print(sequence.sequence)
            ar.append(len(sequence.sequence))
        print(len(sequences))
        print(ar)



    def alikeness(self, s1, s2):
        '''find percentage of similarity between 2 sequences'''
        count = 0
        for i in range(0, len(s1)):
            if s1[i] == s2[i]:
                count += 1
        return count / len(s1)


    def sequencify(self, arr):
        '''wrap string sequences in a Sequence object'''
        for i in range(0, len(arr)):
            arr[i] = Sequence(arr[i])
        return arr


    def cleanSequences(self, seqArray):
        '''clean duplicate sequences in sequence list'''
        returnArr = [seqArray[0].sequence]
        for i in range(1, len(seqArray)):
            if seqArray[i].sequence not in returnArr:
                returnArr.append(seqArray[i].sequence)
        returnArr = self.sequencify(returnArr)
        return returnArr



    def runSimulationSingleAncestor(self, input=None):
        '''
               Simulation starting with a single ancestor.  Will be trying to implement some sort of extinction
               mechanism
        '''
        ancestor = ''
        if len(input.sequence) > 0:
            ancestor = input
        else:
            ancestor = self.getSequences("KT388711")[0]
        print("Simulating")
        #run time in generations
        runtime = 15
        #the rate of reproduction for the influenza A virus is 1.5
        r0 = False
        sequences = [ancestor]
        current = sequences
        newCurrent = []
        while runtime > 0:
            for i in range(0, len(current)):
                newSequence = pm.influenzaMutate(current[i])
                sequences.append(newSequence)
                newCurrent.append(newSequence)
                if r0:
                    newSequence = pm.influenzaMutate(current[i])
                    sequences.append(newSequence)
                    newCurrent.append(newSequence)
            if not r0:
                r0 = True
            else:
                r0 = False
            current = newCurrent
            current = self.cleanSequences(current) #implement or nah?
            newCurrent = []
            runtime -= 1

        #sequences = self.cleanSequences(sequences)
        parseObj.writeToFasta(sequences)
        print("Simulation Complete.")

        self.performMSA()
        self.displayTree()




if __name__ == '__main__':
    pm = ProbabilityModels()
    parseObj = EyeOh()
    seq = input("Please input your own sequence: ")
    simulator = Simulation()
    simulator.runSimulationSingleAncestor(Sequence(seq))











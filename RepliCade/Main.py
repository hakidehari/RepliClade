"""

This is the Main file for Repliclade where the program will be run
All dependencies and classes will be imported here and ran accordingly

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


class Simulation(object):


    def perform_MSA(self):
        '''
                Performs multiple sequence alignment after simulation is run
                using clustalW2, one of the fastest MSA algorithms out there.
                Once finished, displays the sequences and their alignment
        '''
        print("Performing Multiple Sequence Alignment.......")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        clustalw_exe = os.path.join(dir_path, "clustalw2.exe")
        clustalw_cline = ClustalwCommandline(clustalw_exe, infile="influenza.fasta")
        stdout, stderr = clustalw_cline()
        align = AlignIO.read("influenza.aln", "clustal")
        print(align)
        print("Alignment Complete.")


    def display_tree(self):
        '''Creates and displays a phylogenetic tree based on the Multiple Sequence Alignment scores'''
        tree = Phylo.read("influenza.dnd", "newick")
        tree.rooted = True
        Phylo.draw(tree)


    def get_sequences_influenza_B(self, gene):
        '''Gathers sequences from GenBank'''
        con = Connector()
        seq = con.get_gene_data_influenza_B(gene)
        return seq


    def print_sequences(self, sequences):
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


    def clean_sequences(self, seqArray):
        '''clean duplicate sequences in sequence list'''
        returnArr = [seqArray[0].sequence]
        for i in range(1, len(seqArray)):
            if seqArray[i].sequence not in returnArr:
                returnArr.append(seqArray[i].sequence)
        returnArr = self.sequencify(returnArr)
        return returnArr


    def run_influenza_B(self, ancestor, generations):
        '''
               Simulation starting with a single ancestor.  Will be trying to implement some sort of extinction
               mechanism
        '''
        print("Simulating Influenza B genome")
        #run time in generations
        runtime = generations
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
            current = self.clean_sequences(current) #implement or nah?
            newCurrent = []
            runtime -= 1

        #sequences = self.clean_sequences(sequences)
        parseObj.writeToFasta(sequences)
        print("Simulation Complete.")

        self.perform_MSA()
        self.display_tree()


    def run_simulation_input_sequence(self, ancestor, generations, genome_type):
        '''
            Simulation starting with a single ancestor input sequence
        '''
        print("Simulating input genome...")
        


if __name__ == '__main__':
    pm = ProbabilityModels()
    parseObj = EyeOh()
    simulator = Simulation()
    seq = input("Please input your own sequence or leave blank for Influenza B genome: ")
    generations = eval(input("How many generations would you like the program to run for?: "))
    if len(seq) > 0:
        viral_genome_input = input("Is this a viral genome?  Please enter 'Y' for yes and 'N' for no:  ")
        simulator.run_simulation_input_sequence(Sequence(seq), generations, viral_genome_input)
    else:
        influenza_B_segments = simulator.get_sequences_influenza_B(seq)
        influenza_B_segments = [line.replace("\n", "") for line in influenza_B_segments]
        print(influenza_B_segments)
        influenza_B_complete = "".join(influenza_B_segments)
        print(len(influenza_B_complete))
        simulator.run_influenza_B(Sequence(influenza_B_complete), generations)











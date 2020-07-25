"""
This is the main Simulation file.  From here, the simulation will begin
and execute until completion
"""
from PhylogeneticTree import TreeNode
from PhylogeneticTree import PhylogeneticTree
from Generator import Generator
import re
from ProbabilityModels import *

class Simulation(object):

    def __init__(self):
        print("Simulation has started")
        self.virus = ''
        self.rnaught = 0
        self.time_type = ''
        self.is_random_sequence = ''
        self.sequence = ''
        self.generator = Generator()
        self.time_length = 0


    def prompt(self):
        '''
            Prompts the user for information on the sequence to be mutated
        '''

        self.virus = input("Is this genome a virus genome? Enter Y or N:  ")
        while not re.search("[yYnN]", self.virus) or len(self.virus) > 1:
            self.virus = input("Invalid input. Please specify either Y or N for whether the genome sequence is a virus:  ")
        if self.virus.lower() == 'y':
            self.rnaught = input("What is the r-naught value of the virus?:  ")
            try:
                self.rnaught = float(self.rnaught)
            except ValueError:
                while not isinstance(self.rnaught, float):
                    try:
                        self.rnaught = input("Invalid input. Please enter a number for the r-naught value:  ")
                        self.rnaught = float(self.rnaught)
                    except ValueError:
                        continue
        
        self.time_type = input("Would you like the time unit to be in generations or years? Type either G or Y:  ")
        while not re.search("[gGyY]", self.time_type) or len(self.time_type) > 1:
            self.time_type = input("Invalid input. Please enter either G or Y for generations or years:  ")

        self.time_length = input("Please enter a number for the amount of generations or years to run the simulation:  ")
        try:
            self.time_length = int(self.time_length)
        except ValueError:
            while not isinstance(self.time_length, int):
                try:
                    self.time_length = input("Invalid input.  Please enter a valid number for the amount of generations or years to run the simulation:  ")
                    self.time_length = int(self.time_length)
                except ValueError:
                    continue

        self.is_random_sequence = input("Should we generate a randomly generated sequence or will the user provide one? Type either Y or N:  ")
        while not re.search("[yYnN]", self.is_random_sequence) or len(self.is_random_sequence) > 1:
            self.is_random_sequence = input("Invalid input.  Please enter either Y or N if you would like to generate a random sequence:  ")

        if self.is_random_sequence.lower() == 'y':
            seq_len = input("Please enter the desired length for the randomly generated sequence:  ")
            try:
                seq_len = int(seq_len)
            except ValueError:
                while not isinstance(seq_len, int):
                    try:
                        seq_len = input("Invalid input.  Please enter a number for the sequence length:  ")
                        seq_len = int(seq_len)
                    except ValueError:
                        continue
            self.sequence = self.generator.generate_sequence(seq_len)
            #print("The sequence is:\n\n" + self.sequence)


    def simulate(self):
        '''Function which performs the simulation.  Takes the class variables as input taken from the prompts'''
        origin_sequence = TreeNode(sequence=self.sequence)
        nodes = origin_sequence


        for i in range(0, self.time_length):
            if isinstance(nodes, list):
                for k in range (0, len(nodes)):
                    for j in range(0, int(self.rnaught)):
                        new_sequence = TreeNode(sequence=mutate_sequence_random(nodes[k].sequence))
                        nodes[k].add_children(new_sequence)
            else:
                for j in range(0, int(self.rnaught)):
                    new_sequence = TreeNode(sequence=mutate_sequence_random(nodes[k].sequence))
                    nodes.add_children(new_sequence)
            

        #initialize tree
        phylo_tree = PhylogeneticTree(origin_sequence)
        print(len(origin_sequence.children))
        





if __name__ == '__main__':
    sim = Simulation()
    sim.prompt()
    sim.simulate()




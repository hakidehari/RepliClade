
from Sequence import Sequence
import os

class EyeOh:

    def __init__(self):
        pass


    def write_to_fasta(self, seqArray):
        '''Writes sequence to fasta file in Main directory'''
        file = open("influenza.fasta", "w")
        for i in range(0,len(seqArray)):
            file.write(">Sequence%d\n" % float(i+1))
            file.write("%s*" % seqArray[i].sequence)
            file.write("\n\n")
        file.close()


    def check_for_sequence_in_file(self, id):
        '''Check for sequence in file and returns it if exists'''
        try:
            file = open('saved_sequence_data.txt', 'r')
            file.close()
            # Store configuration file values
        except FileNotFoundError:
            file = open("saved_sequence_data.txt", "w")
            file.close()
        existsInFile = False
        file = open("saved_sequence_data.txt", "r")
        for line in file:
            if existsInFile is True:
                return line
            if id in line:
                existsInFile = True
        return None

    def write_sequence_to_file(self, id, seq):
        '''Writes sequence to a text file'''
        file = open("saved_sequence_data.txt", "a")
        file.write(id + "\n" + str(seq) + "\n\n")
        file.close()

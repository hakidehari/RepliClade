
from Sequence import Sequence
import os

class EyeOh:

    def __init__(self):
        pass

    def writeToFasta(self, seqArray):
        file = open("influenza.fasta", "w")
        for i in range(0,len(seqArray)):
            file.write(">Sequence%d\n" % float(i+1))
            file.write("%s*" % seqArray[i].sequence)
            file.write("\n\n")
        file.close()


    def checkForSequenceInFile(self, id):
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

    def writeSequenceToFile(self, id, seq):
        file = open("saved_sequence_data.txt", "a")
        file.write(id + "\n" + str(seq) + "\n\n")
        file.close()

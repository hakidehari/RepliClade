
from Sequence import Sequence

class ParseFasta:

    def __init__(self):
        pass

    def writeToFasta(self, seqArray):
        file = open("influenza.fasta", "w")
        for i in range(0,len(seqArray)):
            file.write(">Sequence%d\n" % i)
            file.write("%s*" % seqArray[i].sequence)
            file.write("\n\n")
        file.close()

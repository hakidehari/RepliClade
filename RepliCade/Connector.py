
from Bio import Entrez
from Bio import SeqIO
from Sequence import Sequence
from EyeOh import EyeOh

class Connector:

    def __init__(self):
        '''initializes Genbank object to connect to Genbank'''
        Entrez.api_key = "0485fcafceed26f42b2d84056582538e9f09"
        Entrez.email = "haki_sed@hotmail.com"
        self.IO = EyeOh()


    def getGeneData(self, gene):
        '''Fetches Gene data from Genbank or from the directory file if it exists in there'''
        print("Fetching influenza sequences..........")
        if gene is None or gene == '':
            f = open("influenzaB.txt", "r")
            segment_ids = f.readlines()
            segments = [segment for segment in segment_ids]
            print(segments)
        else:
            segments = [gene]
            print(segments)
        seqArray = []
        for segment in segments:
            segment_in_file = self.IO.checkForSequenceInFile(segment)
            if segment_in_file is not None:
                seqArray.append(segment_in_file)
                print("Fetched from file")
            else:
                handle = Entrez.efetch(db="nucleotide", id=segment, rettype="gb", retmode="text")
                record = SeqIO.read(handle, "genbank")
                seqArray.append(record.seq)
                handle.close()
                self.IO.writeSequenceToFile(segment, record.seq)
                print("Fetched from NCBI")
        return seqArray


    def printRecord(self):
        print(self.record)


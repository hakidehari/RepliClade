
from Bio import Entrez
from Bio import SeqIO
from Sequence import Sequence
from EyeOh import EyeOh

class Connector:

    def __init__(self):
        Entrez.api_key = "0485fcafceed26f42b2d84056582538e9f09"
        Entrez.email = "haki_sed@hotmail.com"
        self.IO = EyeOh()


    def getGeneData(self, gene):
        print("Fetching influenza sequences..........")
        if gene is None:
            influenza = ["KT388711", "KT388703", "KT388695", "KT388687", "KT388679", "KT388671",
                        "KT388663", "KT388655"]
        else:
            influenza = [gene]
        seqArray = []
        for strain in influenza:
            strainInFile = self.IO.checkForSequenceInFile(strain)
            if strainInFile is not None:
                seqArray.append(Sequence(strainInFile))
                print("Fetched from file")
            else:
                handle = Entrez.efetch(db="nucleotide", id=strain, rettype="gb", retmode="text")
                record = SeqIO.read(handle, "genbank")
                seqArray.append(Sequence(record.seq))
                handle.close()
                self.IO.writeSequenceToFile(strain, record.seq)
                print("Fetched from NCBI")
        return seqArray

    def printRecord(self):
        print(self.record)


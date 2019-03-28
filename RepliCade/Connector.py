
from Bio import Entrez
from Bio import SeqIO
from Sequence import Sequence

class Connector:

    def __init__(self):
        Entrez.api_key = "0485fcafceed26f42b2d84056582538e9f09"
        Entrez.email = "haki_sed@hotmail.com"


    def getGeneData(self):
        print("Fetching influenza sequences..........")
        influenza = ["KT388711", "KT388703", "KT388695", "KT388687", "KT388679", "KT388671",
                     "KT388663", "KT388655"]
        seqArray = []
        for strain in influenza:
            handle = Entrez.efetch(db="nucleotide", id=strain, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            seqArray.append(Sequence(record.seq))
            handle.close()
        print("Fetch Complete.")
        return seqArray

    def printRecord(self):
        print(self.record)


    #KT388711
    #KT388703
    #KT388695
    #KT388687
    #KT388679
    #KT388671
    #KT388663
    #KT388655
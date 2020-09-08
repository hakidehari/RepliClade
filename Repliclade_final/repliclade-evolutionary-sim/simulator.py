from util.genbank_connector import GenBankConnector
from util import file_util
from DNA import genes.GENES


GENE_NAME_LIST = [
    'cytochrome_b'
]

gen_con = GenBankConnector()

class Simulator(object):

    def __init__(self):
        self.gene = ''


    def prompt(self):
        print("Please select a gene from the list so simulate:")
        for gene in GENE_NAME_LIST:
            print(gene)
        gene_choice = input("==>")
        while gene_choice not in GENE_NAME_LIST:
            gene_choice = input("Invalid choice.  Please choose one of the genes listed above: ")
        self.gene = gene_choice

    
    def run_simulation():
        self.prompt()


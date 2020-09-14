from util.genbank_connector import GenBankConnector
from util.file_util import FileStream
from DNA.genes import GENES
from util.seq_util import SequenceUtil
from util.motif_det import DetMotifFinding
import os


GENE_NAME_LIST = [
    'cytochrome_b'
]

gen_con = GenBankConnector()
seq_util = SequenceUtil()
file_util = FileStream()


class Simulator(object):

    def __init__(self):
        self.gene = ''
        self.motif_alg = ''
        self.motif_alg_func = ''


    def gene_prompt(self):
        print("Please select a gene from the list so simulate:")
        for gene in GENE_NAME_LIST:
            print(gene)
        gene_choice = input("==>")
        while gene_choice not in GENE_NAME_LIST:
            gene_choice = input("Invalid choice.  Please choose one of the genes listed above: ")
        self.gene = gene_choice

    
    def prompt_motif_alg(self):
        '''
        Prompt to choose which motif finding algorithm will be used
        on the set of aligned sequences
        '''

        MOTIFS = [
            'deterministic',
            'probabilistic'
        ]

        DET_ALGS = [
            'exhaustive',
            'branch_and_bound',
            'heuristic'
        ]

        PROB_ALGS = [
            ''
        ]

        print('Now that the sequences have been aligned, which off the following motif finding \n \
                    algorithms would you like to implement?')

        for alg in MOTIFS:
            print(alg)
        choice = input("-->")

        while choice.lower() not in MOTIFS:
            print("Invalid input.  Please specify one of the motif finding algorithms provided")
            for alg in MOTIFS:
                print(alg)
            choice = input("-->")

        self.motif_alg = choice

        print("Which type of deterministic algorithm would you like to implement?")

        for alg in DET_ALGS:
            print(alg)
        choice = input("-->")

        while choice.lower() not in DET_ALGS:
            print("Invalid input.  Please specify one of the deterministic motif finding algorithms provided")
            for alg in DET_ALGS:
                print(alg)
            choice = input("-->")

        self.motif_alg_func = choice

    
    def motif_finder_det(self):
        seq_list = file_util.read_from_alignment(self.gene)
        mf_det = DetMotifFinding(seqs=seq_list, size=50)

        if self.motif_alg_func == 'exhaustive':
            sol = mf_det.exhaustive_search()
            print("Exhaustive search solution: ", sol)
        
        if self.motif_alg_func == 'branch_and_bound':
            sol = mf_det.branch_and_bound()
            print("Branch and Bound search solution: ", sol)

        if self.motif_alg_func == 'heuristic':
            sol = mf_det.heuristic_consensus()
            print("Heuristic search solution: ", sol)
    

    def run_simulation(self):
        self.gene_prompt()
        genbank_ids = []
        for gene in GENES[self.gene]:
            genbank_ids.append(gene['id'])
        
        #fetch sequences from genbank
        gene_array = self.fetch_gene_sequence_from_genbank(genbank_ids)

        #write sequences to fasta file
        file_util.write_to_fasta(gene_array, self.gene)

        #align gene sequences
        seq_util.align_sequences_w2(self.gene)

        #prompt motif finding alg
        self.prompt_motif_alg()

        #find motifs within the set of aligned sequences
        if self.motif_alg == 'deterministic':
            self.motif_finder_det()

        if self.motif_alg == 'probabilistic':
            pass

        
    
    def fetch_gene_sequence_from_genbank(self, genes):
        return gen_con.fetch_sequences(genes)


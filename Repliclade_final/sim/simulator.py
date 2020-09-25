from util.genbank_connector import GenBankConnector
from util.file_util import FileStream
from DNA.genes import GENES
from util.seq_util import SequenceUtil
from util.motif_finding import MotifFinding
from util.evolve import JukesCantor, Kimura
import os
import json


GENE_NAME_LIST = [
    'cytochrome_b',
    'cytochrome_b_1'
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


    def file_prompt(self):
        file_prompt = input("Would you like to use an input file?  Please specify Y or N:  ")
        if file_prompt in ['y', 'Y']:
            filename = input("Please specify the name of the file: ")
            return (filename, True)
        return (None, False)

    
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
            'heuristic_stochastic',
            'gibbs'
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

        if self.motif_alg == 'deterministic':
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
        
        if self.motif_alg == 'probabilistic':
            print("Which type of probabilistic algorithm would you like to implement?")

            for alg in PROB_ALGS:
                print(alg)
            choice = input("-->")

            while choice.lower() not in PROB_ALGS:
                print("Invalid input.  Please specift one of the probabilistic motif finding algorithms provided")
                for alg in PROB_ALGS:
                    print(alg)
                choice = input("-->")
            
            self.motif_alg_func = choice


    
    def motif_finder_det(self):
        seq_list = file_util.read_from_alignment()
        mf = MotifFinding(seqs=seq_list, size=50)

        if self.motif_alg_func == 'exhaustive':
            sol = mf.exhaustive_search()
            print("Exhaustive search solution: ", sol)
        
        if self.motif_alg_func == 'branch_and_bound':
            sol = mf.branch_and_bound()
            print("Branch and Bound search solution: ", sol)

        if self.motif_alg_func == 'heuristic':
            sol = mf.heuristic_consensus()
            print("Heuristic search solution: ", sol)

    
    def motif_finder_prb(self):
        seq_list = file_util.read_from_alignment()
        for j in range(8, int(len(seq_list[0])/2)):
            mf = MotifFinding(seqs=seq_list, size=j)

            if self.motif_alg_func == 'heuristic_stochastic':
                sol = mf.heuristic_stochastic()
                print("Heuristic Stochastic search solution: ", sol)
            
            if self.motif_alg_func == 'gibbs':
                sol = mf.gibbs(iterations=100)
                for i in range(len(seq_list)):
                    print(seq_list[i][sol[i]:sol[i]+j])
                print("Gibbs search solution: ", sol)

    
    def prompt_time(self):
        generations = input("Please specify the amount of generations you would like to run the simulation for: ")
        try:
            generations = int(generations)
        except:
            while not isinstance(generations, int):
                generations = input("Invalid input.  Please input the amount of generations you would like to run the simulation for as an integer: ")
                try:
                    generations = int(generations)
                except:
                    continue
        return generations
    

    def simulate(self, c_regions, filename):
        # read sequences from blast
        seq_to_simulate = file_util.read_from_blast_fasta(filename)
        #declare evolutionary model to use
        model = Kimura()
        #have user input the generations
        generations = self.prompt_time()

        generation_dict = {}

        for unit in range(generations):
            current_seqs = []
            for seq in seq_to_simulate:
                new_seq = model.evolve(seq)
                current_seqs.append(new_seq)
            generation_dict[unit] = current_seqs
            seq_to_simulate = current_seqs
        #comment test
        file_util.log_simulation_to_json(generation_dict)
        


    def run_simulation(self):
        file_bool = self.file_prompt()
        if file_bool[1] == False:
            self.gene_prompt()
            genbank_ids = []
            for gene in GENES[self.gene]:
                genbank_ids.append(gene['id'])
            
            #fetch sequences from genbank
            gene_array = self.fetch_gene_sequence_from_genbank(genbank_ids)

            #write sequences to fasta file
            file_util.write_to_fasta(gene_array, self.gene)
            
            #align gene sequences
            seq_util.align_sequences_w2_fasta(self.gene)

            #classify conserved regions of alignment
            conserved_regions = seq_util.calculate_conserved_regions()
            
            print(conserved_regions)
        else:
            filename = file_bool[0]
            #commented out to speed up testing
            #gen_con.run_ncbi_blast_input_file(filename)
            seqs_blast = file_util.read_from_blast(filename)
            file_util.write_to_fasta_blast(seqs_blast, filename)
            seq_util.align_sequences_w2_file(filename)
            conserved_regions = seq_util.calculate_conserved_regions()
            print(filename)

            self.simulate(c_regions=conserved_regions, filename=filename)





        '''
        #prompt motif finding alg
        self.prompt_motif_alg()

        #find motifs within the set of aligned sequences
        if self.motif_alg == 'deterministic':
            self.motif_finder_det()

        if self.motif_alg == 'probabilistic':
            self.motif_finder_prb()
        '''

        
    
    def fetch_gene_sequence_from_genbank(self, genes):
        return gen_con.fetch_sequences(genes)


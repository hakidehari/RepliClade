from util.genbank_connector import GenBankConnector
from util.file_util import FileStream
from DNA.genes import GENES
from util.seq_util import SequenceUtil
from util.motif_finding import MotifFinding
from util.evolve import JukesCantor, Kimura, Felsenstein
import os
import json
import time


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
        print("We have estimated earlier that the total ancestor coalescence time is {}".format(seq_util.coalescence_time))
        print("We can provide an amount of generations for you to run the simulation.  Would you like to enter an amount in years for each generation for this specific sequence?")
        yes_or_no = input("Y or N: ")
        while yes_or_no not in ['N', 'n', 'Y', 'y']:
            yes_or_no = input("Invalid input. Please enter Y or N: ")
        if yes_or_no in ['Y', 'y']:
            years_per_generation = input("Please enter the amount in years for each generation for this input sequence: ")
            try:
                years_per_generation = float(years_per_generation)
            except:
                while not isinstance(years_per_generation, float):
                    years_per_generation = input("Invalid input. Please input a valid number: ")
                    try:
                        years_per_generation = float(years_per_generation)
                    except:
                        continue
            generations = int(seq_util.coalescence_time / years_per_generation)
            print("The simulation will run for {} generations".format(generations))
            return generations
        else:
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

    
    def prompt_model(self):
        evol_models = ["kimura", "jukescantor", "felsenstein"]
        for model in evol_models:
            print(model)
        model = input("Please specify the evolutionary model you would like to use from the ones given above: ")
        while model.lower() not in evol_models:
            model = input("Invalid input.  Please specify the evolutionary model you would like to use from the ones given above: ")
        return model.lower()
        
    
    def cr_prompt(self):
        inp = input("Would you like to consider conserved regions previously identified?  Please input Y or N: ")
        while inp not in ['Y', 'y', 'N', 'n']:
            inp = input("Invalid Input.  Please specify Y or N")
        return inp

    
    def align_results_prompt(self):
        inp = input("Would you like to perform an alignment on the results for the sequence/sequences? Please input 'Y' or 'N': ")
        while inp not in ['Y', 'y', 'N', 'n']:
            inp = input("Invalid Input.  Please specify Y or N")
        return inp

    
    def align_single_or_multiple_prompt(self):
        inp = input("Would you like to align all of the sequences and their results or just a specific sequence? Please specify 'single' or multiple': ")
        while inp.lower() not in ['single', 'multiple']:
            inp = input("Invalid input.  Please specify single or multiple")
        return inp

    
    def which_sequence_align_prompt(self, seq_ids):
        [print(seq_id.lower()) for seq_id in seq_ids]
        seq_ids = [seq.lower() for seq in seq_ids]
        inp = input("Which sequence would you like to align for every iteration throughout the simulation?: ")
        while inp.lower() not in seq_ids:
            inp = input("Invalid input.  Please specify one of the sequence Id's: ")
        return (inp, seq_ids.index(inp.lower()))


    def simulate(self, c_regions, filename):
        # read sequences from blast
        blast_sequences = file_util.read_from_blast_fasta(filename)
        seq_ids = [seq[1] for seq in blast_sequences]
        seq_to_simulate = [seq[0] for seq in blast_sequences]
        
        #have user input the generations
        generations = self.prompt_time()
        evol_model = self.prompt_model()
        cr_inp = self.cr_prompt()

        generation_dict = {}

        #generate simulation model objects
        if evol_model == 'kimura':
            obj_arr = [Kimura() for seq in seq_to_simulate]
        elif evol_model == 'jukescantor':
            obj_arr = [JukesCantor() for seq in seq_to_simulate]
        elif evol_model == 'felsenstein':
            obj_arr = [Felsenstein(seq) for seq in seq_to_simulate]

        #begin simulation
        print("Running Simulation...")
        for unit in range(generations):
            current_seqs = []
            for i in range(len(seq_to_simulate)):
                if cr_inp.lower() == 'y':
                    new_seq = obj_arr[i].evolve_cr(seq_to_simulate[i], c_regions)
                else:
                    new_seq = obj_arr[i].evolve(seq_to_simulate[i])
                current_seqs.append(new_seq)
            generation_dict[unit] = current_seqs
            seq_to_simulate = current_seqs
        print("Simulation Complete.")

        file_util.log_simulation_to_json(generation_dict)
        seq_util.estimate_substitutions_generations(generation_dict, generations)
        if evol_model == "jukescantor":
            seq_util.calculate_divergence_jc(generation_dict, generations)
        if evol_model == "kimura":
            seq_util.calculate_divergence_k2p(generation_dict, generations)
        
        align_results = self.align_results_prompt()

        if align_results.lower() == 'y':
            single_or_multiple = self.align_single_or_multiple_prompt()
            if single_or_multiple == 'single':
                seq_id, index = self.which_sequence_align_prompt(seq_ids)
                seq_util.align_results_w2_single(generation_dict, generations, index, seq_id)
            if single_or_multiple == 'multiple':
                seq_util.align_results_w2_multiple(generation_dict, generations, seq_ids)

    
    def simulate_ancestor(self, sequence, mu, entropy_scores):
        '''
        Simulates using one ancestral sequence inferred
        '''

        generations = self.prompt_time()
        evol_model = self.prompt_model()
        cr_inp = self.cr_prompt()

        print("Beginning simulation...")
        start = time.time()

        generation_dict = {}

        if evol_model == 'kimura':
            model = [Kimura(mu)]
        if evol_model == 'jukescantor':
            model = [JukesCantor(mu)]
        if evol_model == 'felsenstein':
            model = [Felsenstein(sequence)]

        current_seqs = [sequence]
        ext_dict = {}
        dup_dict = {}
        dup_event = False
        ext_event = False
        new_gen = []

        for i in range(generations):
            seq_count = len(current_seqs)
            j = 0
            while j < seq_count:
                dup_event = seq_util.roll_duplication()
                ext_event = seq_util.roll_extinction()
                indel_event = model[j].roll_indel()
                if dup_event:
                    print("Duplication event")
                    new_gen.append(current_seqs[j])
                    new_gen.append(current_seqs[j])
                    dup_dict[i] = "Sequence \n{0}\n was duplicated at time generation {1}".format(current_seqs[j], i)
                    model.append(Kimura(mu) if evol_model == 'kimura' else JukesCantor(mu) if evol_model == 'jukescantor' else Felsenstein(current_seqs[j]) if evol_model == 'felsenstein' else None)
                    dup_event = False
                elif ext_event and seq_count > 1:
                    print("Extinction event")
                    ext_dict[i] = "Sequence \n{0}\n went extinct at time generation {1}".format(current_seqs[j], i)
                    del model[j]
                    del current_seqs[j]
                    j-=1
                    seq_count -= 1
                    ext_event = False
                else:
                    new_seq = model[j].evolve(current_seqs[j])
                    new_gen.append(new_seq)
                j += 1
            #print(len(current_seqs))
            current_seqs = new_gen
            new_gen = []
            
            #break out of simulation if all sequences go extinct
            if len(current_seqs) == 0:
                break
            else:
                generation_dict[i] = current_seqs
        end = time.time()
        print("Simulation Complete.")
        print("Time elapsed: {} seconds".format(end - start))
        print(len(current_seqs))

        if len(current_seqs) == 0:
            print("All sequences went extinct.")
         ######THESE NEED TO BE FIXED #######
        #seq_util.estimate_substitutions_generations(generation_dict, generations)
        #if evol_model == "kimura":
            #seq_util.calculate_divergence_k2p(generation_dict, generations)
        #if evol_model == 'jukescantor':
            #seq_util.calculate_divergence_jc(generation_dict, generations)
        #######################################
        [print(ext_dict[key]) for key in ext_dict]
        [print(dup_dict[key]) for key in dup_dict]
        


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

            print(len(seqs_blast))

            file_util.write_to_fasta_blast(seqs_blast, filename)

            #seq_util.align_sequences_w2_file(filename)

            entropy_scores = seq_util.calculate_conserved_regions()
            
            print(filename)

            aligned_seqs = file_util.read_from_alignment()

            ancestral_seq = seq_util.coalesce(aligned_seqs)

            self.simulate_ancestor(ancestral_seq, seq_util.mu, entropy_scores)





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


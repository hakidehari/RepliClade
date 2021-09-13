from repliclade.util.genbank_connector import GenBankConnector
from repliclade.util.file_util import FileStream
from repliclade.util.seq_util import SequenceUtil
from repliclade.util.evolve import JukesCantor, Kimura, Felsenstein, HKY85
from repliclade.util.evol_tree import TreeNode, print_tree, Phylogenize
from Bio.Phylo.BaseTree import Clade, Tree
import time


gen_con = GenBankConnector()
file_util = FileStream()


class Simulator(object):
    def __init__(self):
        self.seq_util = SequenceUtil()

    @staticmethod
    def file_prompt():
        filename = input("Please specify the name of the file: ")
        return filename

    def prompt_time(self):
        print(
            "We have estimated earlier that the total ancestor coalescence time is {}".format(
                self.seq_util.coalescence_time
            )
        )
        print(
            "We can provide an amount of generations for you to run the simulation.  Would you like to enter an amount in years for each generation for this specific sequence?"
        )
        yes_or_no = input("Y or N: ")
        while yes_or_no not in ["N", "n", "Y", "y"]:
            yes_or_no = input("Invalid input. Please enter Y or N: ")
        if yes_or_no in ["Y", "y"]:
            years_per_generation = input(
                "Please enter the amount in years for each generation for this input sequence: "
            )
            try:
                years_per_generation = float(years_per_generation)
            except:
                while not isinstance(years_per_generation, float):
                    years_per_generation = input(
                        "Invalid input. Please input a valid number: "
                    )
                    try:
                        years_per_generation = float(years_per_generation)
                    except:
                        continue
            generations = int(self.seq_util.coalescence_time / years_per_generation)
            print("The simulation will run for {} generations".format(generations))
            return generations
        else:
            generations = input(
                "Please specify the amount of generations you would like to run the simulation for: "
            )
            try:
                generations = int(generations)
            except:
                while not isinstance(generations, int):
                    generations = input(
                        "Invalid input.  Please input the amount of generations you would like to run the simulation for as an integer: "
                    )
                    try:
                        generations = int(generations)
                    except:
                        continue
            return generations

    @staticmethod
    def prompt_model():
        evol_models = ["kimura", "jukescantor", "felsenstein", "hasegawa"]
        for model in evol_models:
            print(model)
        model = input(
            "Please specify the evolutionary model you would like to use from the ones given above: "
        )
        while model.lower() not in evol_models:
            model = input(
                "Invalid input.  Please specify the evolutionary model you would like to use from the ones given above: "
            )
        return model.lower()

    @staticmethod
    def cr_prompt():
        inp = input(
            "Would you like to consider conserved regions previously identified?  Please input Y or N: "
        )
        while inp not in ["Y", "y", "N", "n"]:
            inp = input("Invalid Input.  Please specify Y or N")
        return inp

    @staticmethod
    def align_results_prompt():
        inp = input(
            "Would you like to perform an alignment on the results for the sequence/sequences? Please input 'Y' or 'N': "
        )
        while inp not in ["Y", "y", "N", "n"]:
            inp = input("Invalid Input.  Please specify Y or N")
        return inp

    @staticmethod
    def align_single_or_multiple_prompt():
        inp = input(
            "Would you like to align all of the sequences and their results or just a specific sequence? Please specify 'single' or multiple': "
        )
        while inp.lower() not in ["single", "multiple"]:
            inp = input("Invalid input.  Please specify single or multiple")
        return inp

    @staticmethod
    def which_sequence_align_prompt(seq_ids):
        [print(seq_id.lower()) for seq_id in seq_ids]
        seq_ids = [seq.lower() for seq in seq_ids]
        inp = input(
            "Which sequence would you like to align for every iteration throughout the simulation?: "
        )
        while inp.lower() not in seq_ids:
            inp = input("Invalid input.  Please specify one of the sequence Id's: ")
        return (inp, seq_ids.index(inp.lower()))

    def simulate_ancestor(self, sequence, mu, entropy_scores, model=None):
        """
        Simulates using one ancestral sequence inferred
        """

        generations = 50000
        evol_model = model
        cr_inp = "N"
        # generations = self.prompt_time()
        # evol_model = self.prompt_model()
        # cr_inp = self.cr_prompt()

        print("Beginning simulation...")
        start = time.time()

        generation_dict = {}

        if evol_model == "kimura":
            model = [Kimura(mu)]
        if evol_model == "jukescantor":
            model = [JukesCantor(mu)]
        if evol_model == "felsenstein":
            model = [Felsenstein(sequence)]
        if evol_model == "hasegawa":
            model = [HKY85(sequence)]

        node_num = 1
        clades = [Clade(name=str(node_num))]
        tree = [TreeNode(sequence, root=True, num=node_num)]
        current_seqs = [sequence]
        ext_dict = {}
        dup_dict = {}
        dup_event = False
        ext_event = False
        new_gen = []
        all_seqs = []
        all_nodes = []
        all_clades = []

        for i in range(generations):
            seq_count = len(current_seqs)
            j = 0
            while j < seq_count:
                dup_event = self.seq_util.roll_duplication()
                ext_event = self.seq_util.roll_extinction()
                indel_event = self.seq_util.roll_indel()
                current_seqs[j] = model[j].evolve(current_seqs[j])
                tree[j].set_sequence(current_seqs[j])
                if dup_event:
                    print(f"Duplication event at time generation {i}")
                    node_num += 1
                    new_gen.append(current_seqs[j])
                    new_gen.append(current_seqs[j])
                    new_clade = Clade(name=str(node_num))
                    clades.append(new_clade)
                    clades[j].clades.append(new_clade)
                    child = TreeNode(current_seqs[j], num=node_num)
                    tree.append(child)
                    tree[j].add_child(child)
                    dup_dict[
                        i
                    ] = "Sequence \n{0}\n was duplicated at time generation {1}".format(
                        current_seqs[j], i
                    )
                    model.append(
                        Kimura(mu)
                        if evol_model == "kimura"
                        else JukesCantor(mu)
                        if evol_model == "jukescantor"
                        else Felsenstein(current_seqs[j])
                        if evol_model == "felsenstein"
                        else HKY85(current_seqs[j])
                        if evol_model == "hasegawa"
                        else None
                    )
                    dup_event = False
                elif ext_event and seq_count > 1:
                    print(f"Extinction event at time generation {i}")
                    ext_dict[
                        i
                    ] = "Sequence \n{0}\n went extinct at time generation {1}".format(
                        current_seqs[j], i
                    )
                    all_seqs.append(current_seqs[j])
                    all_nodes.append(tree[j])
                    all_clades.append(clades[j])
                    del model[j]
                    del current_seqs[j]
                    del tree[j]
                    del clades[j]
                    j -= 1
                    seq_count -= 1
                    ext_event = False
                # elif indel_event:
                # current_seqs[j] = model[j].execute_indel(current_seqs[j])
                # tree[j].set_sequence(current_seqs[j])
                # new_gen.append(current_seqs[j])
                else:
                    new_gen.append(current_seqs[j])

                j += 1
            # print(len(current_seqs))
            current_seqs = new_gen
            new_gen = []

            # break out of simulation if all sequences go extinct
            if len(current_seqs) == 0:
                break
            else:
                generation_dict[i] = current_seqs

        all_clades.extend(clades)
        all_seqs.extend(current_seqs)
        all_nodes.extend(tree)
        end = time.time()
        print("Simulation Complete.")
        print("Time elapsed: {} seconds".format(end - start))
        print(len(current_seqs))

        if len(current_seqs) == 0:
            print("All sequences went extinct.")

        [print(ext_dict[key]) for key in ext_dict]
        [print(dup_dict[key]) for key in dup_dict]

        full_output = {**ext_dict, **dup_dict}

        file_util.log_simulation_output_to_json(full_output)

        all_nodes.sort(key=lambda x: x.num)

        for node in all_nodes:
            print("Here is the tree for node {0}:\n".format(node.num))
            print_tree(node)

        root_clade = [clade for clade in all_clades if clade.name == "1"][0]
        true_tree = Tree(root=root_clade)
        print(true_tree)
        # Phylo.draw(true_tree)

        return all_nodes, all_seqs, true_tree

    def run_simulation(self):
        filename = self.file_prompt()

        # commented out to speed up testing
        gen_con.run_ncbi_blast_input_file(filename)

        seqs_blast = file_util.read_from_blast(filename)

        print(len(seqs_blast))

        file_util.write_to_fasta_blast(seqs_blast, filename)

        self.seq_util.align_sequences_muscle_file(filename)

        entropy_scores = self.seq_util.calculate_conserved_regions()

        print(filename)

        average_entropy = sum(entropy_scores[key] for key in entropy_scores) / len(
            entropy_scores
        )
        print(f"Average entropy over all of the sequences: {average_entropy}")

        aligned_seqs = file_util.read_from_alignment()

        ancestral_seq = self.seq_util.coalesce_v2(aligned_seqs)

        tree, post_sim_seqs = self.simulate_ancestor(
            ancestral_seq, self.seq_util.mu, entropy_scores
        )

        filename_results = file_util.write_to_fasta_sim_results(
            post_sim_seqs, tree, filename
        )

        self.seq_util.align_sequences_muscle_file(filename_results)

        sim_aligned_seqs = file_util.read_from_alignment_results()

        # print(sim_aligned_seqs)

        # sim_aligned_seqs.sort(key=lambda x: int(x[1]))

        print([seq[1] for seq in sim_aligned_seqs])

        # print(sim_aligned_seqs)

        phylo = Phylogenize([seq[0] for seq in sim_aligned_seqs])

        # print(phylo.calculate_distance_k2p())
        tree_prompt = phylo.prompt_tree_builder()

        if tree_prompt.lower() == "upgma" or tree_prompt == "nj":
            calculator, dm = phylo.biopython_calc_distances_upgma_nj()
            phylo.build_tree_upgma_nj(calculator, dm, tree_prompt)
        if tree_prompt.lower() == "parsimony":
            phylo.build_tree_parsimony()
        if tree_prompt.lower() == "ml":
            phylo.maximum_likelihood()

    @staticmethod
    def fetch_gene_sequence_from_genbank(genes):
        return gen_con.fetch_sequences(genes)

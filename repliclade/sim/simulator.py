from repliclade.util.genbank_connector import GenBankConnector
from repliclade.util.file_util import FileStream
from repliclade.util.seq_util import SequenceUtil
from repliclade.util.evolve import JukesCantor, Kimura, Felsenstein, HKY85
from repliclade.util.evol_tree import TreeNode, print_tree, Phylogenize
from repliclade.util.prompts import prompt_model, file_prompt, prompt_time
from Bio.Phylo.BaseTree import Clade, Tree
import time


GEN_CON = GenBankConnector()
FILE_UTIL = FileStream()
SEQ_UTIL = SequenceUtil()


class Simulator(object):
    def __init__(self):
        self.mutable_model_mapping = {
            JukesCantor: 1,
            Kimura: 1,
            Felsenstein: 1,
            HKY85: 1,
        }

    def simulate_ancestor(self, sequence, mu, entropy_scores, model=None):
        """
        Simulates using one ancestral sequence inferred
        """

        self.mutable_model_mapping.update(
            {JukesCantor: mu, Kimura: mu, Felsenstein: sequence, HKY85: sequence}
        )
        
        generations = prompt_time(SEQ_UTIL.coalescence_time)
        evol_model = prompt_model()
        # cr_inp = self.cr_prompt() TODO

        print("Beginning simulation...")
        start = time.time()

        generation_dict = dict()
        sim_seqs = list()

        if evol_model == "kimura":
            model = Kimura
        if evol_model == "jukescantor":
            model = JukesCantor
        if evol_model == "felsenstein":
            model = Felsenstein
        else:
            model = HKY85

        sim_seqs.append(model(self.mutable_model_mapping[model]))

        node_num = 1
        clades = [Clade(name=str(node_num))]
        tree = [TreeNode(sequence, root=True, num=node_num)]
        current_seqs = [sequence]
        ext_dict = dict()
        dup_dict = dict()
        dup_event = False
        ext_event = False
        new_gen = list()
        all_seqs = list()
        all_nodes = list()
        all_clades = list()

        for i in range(generations):
            seq_count = len(current_seqs)
            j = 0
            while j < seq_count:
                dup_event = SEQ_UTIL.roll_duplication()
                ext_event = SEQ_UTIL.roll_extinction()
                # indel_event = SEQ_UTIL.roll_indel() TODO
                current_seqs[j] = sim_seqs[j].evolve(current_seqs[j])
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

                    if model in {Felsenstein, HKY85}:
                        self.mutable_model_mapping.update({model: current_seqs[j]})

                    sim_seqs.append(model(self.mutable_model_mapping[model]))
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
                    del sim_seqs[j]
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

        FILE_UTIL.log_simulation_output_to_json(full_output)

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
        filename = file_prompt()

        # commented out to speed up testing
        # gen_con.run_ncbi_blast_input_file(filename)

        seqs_blast = FILE_UTIL.read_from_blast(filename)

        FILE_UTIL.write_to_fasta_blast(seqs_blast, filename)

        SEQ_UTIL.align_sequences_muscle_file(filename)

        entropy_scores = SEQ_UTIL.calculate_conserved_regions()

        print(filename)

        average_entropy = sum(entropy_scores[key] for key in entropy_scores) / len(
            entropy_scores
        )
        print(f"Average entropy over all of the sequences: {average_entropy}")

        aligned_seqs = FILE_UTIL.read_from_alignment()

        ancestral_seq = SEQ_UTIL.coalesce_v2(aligned_seqs)

        all_nodes, post_sim_seqs, tree = self.simulate_ancestor(
            ancestral_seq, SEQ_UTIL.mu, entropy_scores
        )

        filename_results = FILE_UTIL.write_to_fasta_sim_results(
            post_sim_seqs, all_nodes, filename
        )

        SEQ_UTIL.align_sequences_muscle_file(filename_results)

        sim_aligned_seqs = FILE_UTIL.read_from_alignment_results()

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
        return GEN_CON.fetch_sequences(genes)

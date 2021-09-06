from __future__ import division
from Bio.Align.Applications import (
    ClustalwCommandline,
    MuscleCommandline,
)
from Bio import Phylo, AlignIO
from Bio.Seq import Seq
from datetime import datetime
from repliclade.util.file_util import FileStream
import numpy as np
import os
import random
import time

file_tool = FileStream()


class SequenceUtil(object):
    """
    All functions related to working directly with sequences other than the motif
    finding functions which are contained in their own files
    """

    gens_passed_dup = 0
    gens_passed_ext = 0
    gens_passed_indel = 0

    def align_sequences_w2_fasta(self, gene_name):
        """
        Aligns sequences using the ClustalW2 executable

        input: gene_name --> to specify which fasta file to read from
        output:  Displays a phylogenetic tree
        """

        print("Aligning sequences using ClustalW2...")
        time = self.get_time()
        in_file = (
            os.getcwd()
            + os.path.sep
            + "alignment"
            + os.path.sep
            + "fastas"
            + os.path.sep
            + "{0}.fasta".format(gene_name)
        )
        out_file = (
            os.getcwd()
            + os.path.sep
            + "alignment"
            + os.path.sep
            + "align"
            + os.path.sep
            + "{0}_out_{1}.aln".format(gene_name, time)
        )

        if os.name == "nt":
            executable = (
                os.getcwd()
                + os.path.sep
                + "alignment"
                + os.path.sep
                + "executables"
                + os.path.sep
                + "clustalw2.exe"
            )
        else:
            executable = (
                os.getcwd()
                + os.path.sep
                + "alignment"
                + os.path.sep
                + "executables"
                + os.path.sep
                + "clustalw2"
            )

        clustalw_cline = ClustalwCommandline(
            executable, infile=in_file, outfile=out_file
        )
        stdout, stderr = clustalw_cline()
        align = AlignIO.read(out_file, "clustal")

        print("Finished aligning sequences.")
        print(align)

        self.display_phylo_tree(gene_name, in_file)

    def align_sequences_w2_file(self, filename):
        print("Aligning sequences using ClustalW2...")
        print("file name:", filename)
        time = self.get_time()
        in_file = (
            os.getcwd()
            + os.path.sep
            + "DNA"
            + os.path.sep
            + "{}.fasta".format(filename)
        )
        out_file = (
            os.getcwd()
            + os.path.sep
            + "alignment"
            + os.path.sep
            + "align"
            + os.path.sep
            + "{0}_out_{1}.aln".format(filename, time)
        )

        if os.name == "nt":
            executable = (
                os.getcwd()
                + os.path.sep
                + "alignment"
                + os.path.sep
                + "executables"
                + os.path.sep
                + "clustalw2.exe"
            )
        else:
            executable = (
                os.getcwd()
                + os.path.sep
                + "alignment"
                + os.path.sep
                + "executables"
                + os.path.sep
                + "clustalw2"
            )

        clustalw_cline = ClustalwCommandline(
            executable, infile=in_file, outfile=out_file
        )
        stdout, stderr = clustalw_cline()
        align = AlignIO.read(out_file, "clustal")

        print("Finished aligning sequences.")
        print(align)

        self.display_phylo_tree(filename)

    def align_sequences_muscle_file(self, filename):
        print("Aligning sequences using Muscle...")
        print("file name:", filename)
        time = self.get_time()
        in_file = (
            os.getcwd()
            + os.path.sep
            + "DNA"
            + os.path.sep
            + "{}.fasta".format(filename)
        )
        out_file = (
            os.getcwd()
            + os.path.sep
            + "alignment"
            + os.path.sep
            + "align"
            + os.path.sep
            + "{0}_out_{1}.aln".format(filename, time)
        )

        if os.name == "nt":
            executable = (
                os.getcwd()
                + os.path.sep
                + "alignment"
                + os.path.sep
                + "executables"
                + os.path.sep
                + "muscle.exe"
            )
        else:
            executable = (
                os.getcwd()
                + os.path.sep
                + "alignment"
                + os.path.sep
                + "executables"
                + os.path.sep
                + "muscle"
            )

        muscle_cline = MuscleCommandline(
            executable, input=in_file, out=out_file, clw=True
        )
        stdout, stderr = muscle_cline()
        align = AlignIO.read(out_file, "clustal")

        print("Finished aligning sequences.")
        print(align)

        self.display_phylo_tree(align)

    def align_results_w2_single(self, generation_dict, generations, index, seq_id):
        print("Aligning results...")
        sequences_gens = [generation_dict[i][index] for i in range(generations)]
        file_tool.write_to_fasta_results(seq_id, sequences_gens)
        time = self.get_time()
        in_file = file_tool.most_recent_fasta_results()
        out_file = (
            os.getcwd()
            + os.path.sep
            + "alignment"
            + os.path.sep
            + "align"
            + os.path.sep
            + "{0}_{1}_aln".format(seq_id.replace("|", ""), time)
        )

        if os.name == "nt":
            executable = (
                os.getcwd()
                + os.path.sep
                + "alignment"
                + os.path.sep
                + "executables"
                + os.path.sep
                + "clustalw2.exe"
            )
        else:
            executable = (
                os.getcwd()
                + os.path.sep
                + "alignment"
                + os.path.sep
                + "executables"
                + os.path.sep
                + "clustalw2"
            )

        clustalw_cline = ClustalwCommandline(
            executable, infile=in_file, outfile=out_file
        )
        stdout, stderr = clustalw_cline()
        align = AlignIO.read(out_file, "clustal")

        print("Finished aligning sequences.")
        print(align)

        self.display_phylo_tree_results()

    def align_results_w2_multiple(self, generation_dict, generations, seq_ids):
        print("Aligning results...")
        sequences = generation_dict[generations - 1]
        file_tool.write_to_fasta_results_multiple(seq_ids, sequences)
        time = self.get_time()
        in_file = file_tool.most_recent_fasta_results()
        out_file = (
            os.getcwd()
            + os.path.sep
            + "alignment"
            + os.path.sep
            + "align"
            + os.path.sep
            + "post_sim_alignment_{}.aln".format(time)
        )

        if os.name == "nt":
            executable = (
                os.getcwd()
                + os.path.sep
                + "alignment"
                + os.path.sep
                + "executables"
                + os.path.sep
                + "clustalw2.exe"
            )
        else:
            executable = (
                os.getcwd()
                + os.path.sep
                + "alignment"
                + os.path.sep
                + "executables"
                + os.path.sep
                + "clustalw2"
            )

        clustalw_cline = ClustalwCommandline(
            executable, infile=in_file, outfile=out_file
        )
        stdout, stderr = clustalw_cline()
        align = AlignIO.read(out_file, "clustal")

        print("Finished aligning sequences.")
        print(align)

        self.display_phylo_tree_results()

    def display_phylo_tree(self, alignment):
        """
        Reads from the dnd file produced by the ClustalW2 alignment
        and displays the phylogenetic tree based on the results of the alignment

        input: gene_name --> to specify the alignment file
               filename  --> specifies the path to the file
        output: Diplays phylogenetic tree using the Phylo biopython module
        """
        from util.evol_tree import Phylogenize

        # MP
        print("Displaying Phylogenetic tree of sequences...")
        """scorer = ParsimonyScorer()
        searcher = NNITreeSearcher(scorer)
        constructor = ParsimonyTreeConstructor(searcher)
        pars_tree = constructor.build_tree(alignment)"""

        # UPGMA
        """phylo = Phylogenize()
        calc, dm = phylo.biopython_calc_distances_upgma_nj()
        phylo.build_tree_upgma_nj(calc, dm, 'upgma')"""

        # NJ
        phylo = Phylogenize()
        calc, dm = phylo.biopython_calc_distances_upgma_nj()
        phylo.build_tree_upgma_nj(calc, dm, "nj")

    def display_phylo_tree_results(self):
        print("Displaying Phylogenetic tree of sequences...")
        dnd_file = file_tool.most_recent_dnd_results()
        tree = Phylo.read(dnd_file, "newick")
        Phylo.draw(tree)

    def calculate_conserved_regions(self):
        """
        Naively determines conserved regions based on the alignment of the sequences

        input: no input but takes the most recent alignment file
        output: shannons entropy score for each column stored in a dictionary
        """
        # fetch aligned sequences from most recent alignment
        aligned_seqs = file_tool.read_from_alignment()

        # shannons entropy score for each col
        shannon_entropy_dict = {}

        seq_len = len(aligned_seqs[0])
        total_seqs = len(aligned_seqs)
        frequency_dict = {}
        potential_values = ["A", "G", "T", "C", "-"]

        # iterate over sequence length
        for i in range(seq_len):
            # initialize nuc count dict
            nuc_count = {}
            # iterate over sequences
            for seq in aligned_seqs:
                if seq[i] in nuc_count:
                    nuc_count[seq[i]] += 1
                else:
                    nuc_count[seq[i]] = 1
            # set default values for any nucs not counted to 1*10e-6
            for nuc in potential_values:
                if nuc not in nuc_count:
                    nuc_count[nuc] = 0.000001

            shannon_entropy_dict[i] = -1 * sum(
                (nuc_count[key] / total_seqs) * np.log2(nuc_count[key] / total_seqs)
                for key in nuc_count
            )

        return shannon_entropy_dict

    def check_for_cr(self, index, c_regions):
        """
        checks if the nucleotide is a conserved region.  If it is, does not mutate
        """
        for region in c_regions:
            start_index = c_regions[region][0]
            end_index = c_regions[region][1]
            if index >= start_index and index <= end_index:
                return True
        return False

    def all_same(self, items):
        return all(x == items[0] for x in items)

    def sequencify(self, seq_array):
        return [Seq(seq) for seq in seq_array]

    def get_time(self):
        return datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

    def roll_duplication(self):
        alpha = 0.00000001
        prb = alpha * self.gens_passed_dup
        roll = random.random()
        if roll < prb:
            self.gens_passed_dup = 0
            return True
        self.gens_passed_dup += 1
        return False

    def roll_extinction(self):
        alpha = 0.00000001
        prb = alpha * self.gens_passed_ext
        roll = random.random()
        if roll < prb:
            self.gens_passed_ext = 0
            return True
        self.gens_passed_ext += 1
        return False

    def roll_indel(self):
        # Julienne M. Mullaney1,2, Ryan E. Mills4, W. Stephen Pittard5 and Scott E. Devine1,2,3,∗
        THRESHOLD = 0.2
        roll = random.random()
        if roll < THRESHOLD:
            return True
        return False

    def get_gc_content(self, seq):
        if len(seq) == 0:
            return 0
        return float((seq.count("G") + seq.count("C")) / len(seq))

    def get_nuc_count(self, seq):
        count_dict = {"A": 0, "G": 0, "T": 0, "C": 0}
        for key in count_dict:
            count_dict[key] = seq.count(key)
        return count_dict

    def get_nuc_count_multiple(self, seqs):
        count_dict = {"A": 0, "G": 0, "T": 0, "C": 0}
        for seq in seqs:
            for key in count_dict:
                count_dict[key] += seq.count(key)
            return count_dict

    def prompt_theta_method(self):
        theta_method = input(
            "Please enter a choice of method for estimation of θ (variation)\nYou can choose either the Watterson Method (watterson)(Watterson 1975) or the Fay and Wu Method (faywu)(Fay and Wu 2000) "
        )
        while theta_method.lower() not in ["watterson", "faywu"]:
            theta_method = input(
                "Invalid input.  Please enter either watterson or faywu: "
            )
        return theta_method

    def estimate_substitutions_generations(self, generation_dict, generations):
        first_gen = generation_dict[0]
        final_gen = generation_dict[generations - 1]

        # calculate divergence
        for i in range(len(first_gen)):
            score = 0
            origin_seq = first_gen[i]
            final_seq = final_gen[i]
            for j in range(len(origin_seq)):
                if origin_seq[j] == final_seq[j]:
                    score += 1
            print(
                "Sequence {0} in the first and last generation are {1} percent similar".format(
                    i, float(score / len(origin_seq)) * 100
                )
            )

    def coalesce_v2(self, sequences):

        from util.evolve import JukesCantor, Kimura, Felsenstein, HKY85
        from sim.simulator import Simulator
        import math

        def merge_sequences(seq1, seq2):
            ancestor = ""
            for i in range(len(seq1)):
                s1 = seq1[i]
                s2 = seq2[i]

                if s1 == s2 and s1 != "-":
                    ancestor += s1
                elif s1 == "-" and s2 == "-":
                    ancestor += random.choice(["A", "G", "T", "C"])
                elif s1 == "-" or s2 == "-":
                    ancestor += s1 if s1 != "-" else s2
                else:
                    nuc = random.choice([s1, s2])
                    ancestor += nuc
            return ancestor

        total_seqs = len(sequences)

        theta_method = self.prompt_theta_method()

        if theta_method == "watterson":
            print("Estimating effective population size using the Watterson method...")

            eff_pop_size = self.estimate_eff_pop_size_watterson(sequences)

            time_to_coalescence = sum(
                (4 * eff_pop_size) / (i * (i - 1)) for i in range(2, total_seqs + 1)
            )
        else:
            print("Estimating effective population size using the Fay and Wu method...")

            eff_pop_size = self.estimate_eff_pop_size_faywu(sequences)

            time_to_coalescence = sum(
                (4 * eff_pop_size) / (i * (i - 1)) for i in range(2, total_seqs + 1)
            )

        sim_obj = Simulator()

        print("Please select an evolutionary model for the coalescent simulation.\n")

        evol_model = sim_obj.prompt_model()

        print("Simulating the Coalescent.  This may take a while...")

        start = time.time()

        seq_len = len(sequences)

        total_generations = 0

        while seq_len > 1:
            rand1 = random.randint(0, len(sequences) - 1)
            rand2 = random.randint(0, len(sequences) - 1)

            while rand1 == rand2:
                rand1 = random.randint(0, len(sequences) - 1)
                rand2 = random.randint(0, len(sequences) - 1)

            seq1 = sequences[rand1]
            seq2 = sequences[rand2]

            temp1 = rand1 if rand1 < rand2 else rand2
            temp2 = rand2 if rand2 > rand1 else rand1

            sequences = (
                sequences[:temp1]
                + sequences[temp1 + 1 : temp2]
                + sequences[temp2 + 1 :]
            )

            merge = False
            ancestor = ""
            i = 1.0

            scaled_factor = 100.0

            nc2 = float((seq_len * (seq_len - 1)) / 2)

            model1 = (
                JukesCantor(self.mu)
                if evol_model == "jukescantor"
                else Kimura(self.mu)
                if evol_model == "kimura"
                else Felsenstein(seq1)
                if evol_model == "felsenstein"
                else HKY85(seq1)
                if evol_model == "hasegawa"
                else None
            )
            model2 = (
                JukesCantor(self.mu)
                if evol_model == "jukescantor"
                else Kimura(self.mu)
                if evol_model == "kimura"
                else Felsenstein(seq2)
                if evol_model == "felsenstein"
                else HKY85(seq2)
                if evol_model == "hasegawa"
                else None
            )

            while not merge:
                prob_coalesce = (
                    ((1.0 - (1.0 / (nc2 * (2.0 * eff_pop_size)))) ** (i - 1))
                    * (nc2 * (1.0 / (2.0 * eff_pop_size)))
                    * scaled_factor
                )  # scaled factor is correction factor to speed up coalescent simulation
                # prob_coalesce = (1/(2*eff_pop_size)) * math.e**((-1*(i-1))/(2*eff_pop_size))

                coal_rand = random.random()

                if coal_rand > prob_coalesce:
                    seq1 = model1.evolve(seq1)
                    seq2 = model2.evolve(seq2)
                # maximum value
                else:
                    ancestor = merge_sequences(seq1, seq2)
                    merge = True
                    print(
                        "Common Ancestor found {} generations back.".format(
                            i * scaled_factor
                        )
                    )
                i += 1.0
            total_generations += (i - 1) * scaled_factor
            sequences.append(ancestor)
            # eff_pop_size = self.estimate_eff_pop_size_watterson_no_input_mu(sequences)
            seq_len -= 1
        end = time.time()

        print(
            "Coalescent Simulation Complete.  The simulation took {} seconds".format(
                end - start
            )
        )
        print("Ancestral Sequence Inferred:\n", sequences[0])
        print(
            "Took approximately {} generations for full coalescence of all sequences".format(
                total_generations
            )
        )
        return sequences[0]

    def coalesce(self, sequences):
        """
        Infers the ancestral sequence of a group of sequences generated from BLAST
        input:  array of DNA sequences
        output:  inferred ancestral sequence
        """
        total_seqs = len(sequences)
        sequences_clone = sequences

        def merge(seq1, seq2):
            ancestor = ""
            for i in range(len(seq1)):
                s1 = seq1[i]
                s2 = seq2[i]

                if s1 == s2 and s1 != "-":
                    ancestor += s1
                elif s1 == "-" and s2 == "-":
                    ancestor += random.choice(["A", "G", "T", "C"])
                elif s1 == "-" or s2 == "-":
                    ancestor += s1 if s1 != "-" else s2
                else:
                    nuc = random.choice([s1, s2])
                    ancestor += nuc
            return ancestor

        while len(sequences) > 1:

            rand1 = random.randint(0, len(sequences) - 1)
            rand2 = random.randint(0, len(sequences) - 1)

            while rand1 == rand2:
                rand1 = random.randint(0, len(sequences) - 1)
                rand2 = random.randint(0, len(sequences) - 1)

            seq1 = sequences[rand1]
            seq2 = sequences[rand2]

            temp1 = rand1 if rand1 < rand2 else rand2
            temp2 = rand2 if rand2 > rand1 else rand1

            sequences = (
                sequences[:temp1]
                + sequences[temp1 + 1 : temp2]
                + sequences[temp2 + 1 :]
            )

            new_seq = merge(seq1, seq2)

            sequences.append(new_seq)

        print("Ancestral sequence inferred: ", sequences[0])

        theta_method = self.prompt_theta_method()

        if theta_method == "watterson":
            print("Estimating effective population size using the Watterson method...")

            eff_pop_size = self.estimate_eff_pop_size_watterson(sequences_clone)

            time_to_coalescence = sum(
                (4 * eff_pop_size) / (i * (i - 1)) for i in range(2, total_seqs + 1)
            )
        else:
            print("Estimating effective population size using the Fay and Wu method...")

            eff_pop_size = self.estimate_eff_pop_size_faywu(sequences_clone)

            time_to_coalescence = sum(
                (4 * eff_pop_size) / (i * (i - 1)) for i in range(2, total_seqs + 1)
            )

        print(
            "These sequences shared a common ancestor roughly {} years ago.".format(
                time_to_coalescence
            )
        )
        return sequences[0]

    def estimate_eff_pop_size_using_N(self):
        """
        Estimates effective population size
        """
        census = ""
        try:
            census = int(
                input(
                    "Please enter a number for the census population size estimate:  "
                )
            )
        except:
            while not isinstance(census, int):
                census = input(
                    "invalid input.  Please enter a number for the census population size estimate:  "
                )
                try:
                    census = int(census)
                except:
                    continue
        nf = int(census / 2)
        nm = int(census / 2)
        return (4 * nm * nf) / (nm + nf)

    def estimate_eff_pop_size_jf_phylo(self, sequences):
        """
        Estimates the effective population size using Joseph Felsenstein's method
        of phylogenetic estimates
        """
        pass

    def promt_mutation_rate(self):
        """
        Prompts the user for mutation rate
        """
        inp = input(
            "Would you like to estimate the effective pop size with our own input? Please enter y or n: "
        )
        while inp.lower() not in ["y", "n"]:
            inp = input("Invalid input.  Please enter y or n: ")
        if inp.lower() == "y":
            mu = 0
            try:
                mu = float(input("Please enter a value for the mutation rate: "))
            except:
                while not isinstance(mu, float):
                    mu = input(
                        "Invalid input.  Please enter a valid number for the mutation rate"
                    )
                    try:
                        mu = float(mu)
                    except:
                        continue
            return mu
        else:
            return None

    def estimate_eff_pop_size_faywu(self, sequences):
        """
        Estimates the effective population size using the Fay and Wu estimator

        Input: array of sequences
        Output: effective population size
        """

        total_seqs = len(sequences)
        seq_len = len(sequences[0])
        # Mu - most often has the rate of 10e-4.  So we will use that value here
        mu = 0.00001
        # number of segregating sites
        K = 0

        threshold = 0.9
        tracking_dict = {}
        final_dict = {}

        for k in range(total_seqs):
            for i in range(seq_len):
                for seq in sequences[: k + 1]:
                    if seq[i] in tracking_dict:
                        tracking_dict[seq[i]] += 1
                    else:
                        tracking_dict[seq[i]] = 1
                highest_char_cnt = -1
                for key in tracking_dict:
                    if tracking_dict[key] > highest_char_cnt:
                        highest_char_cnt = tracking_dict[key]
                if highest_char_cnt / (k + 1) < threshold:
                    K += 1
                tracking_dict = {}
            final_dict[k + 1] = K
            K = 0

        theta_h = sum(
            ((i ** 2) * final_dict[i]) / ((total_seqs * (total_seqs - 1)) / 2)
            for i in range(1, total_seqs)
        )

        print("theta_h: ", theta_h)
        print("Default μ: ", mu)

        self.mu = mu
        Ne = theta_h / (4 * mu)
        print(
            "Effective Population size using Fay and Wu estimator with default μ: ", Ne
        )

        time_to_coalescence = sum(
            (4 * Ne) / (i * (i - 1)) for i in range(2, total_seqs + 1)
        )
        self.coalescence_time = time_to_coalescence
        print(
            "Coalescence time using Effective Population size from default μ: ",
            time_to_coalescence,
        )

        inp_mu = self.promt_mutation_rate()

        if inp_mu is not None:
            Ne = theta_h / (4 * inp_mu)
            print(
                "Effective Population size using Fay and Wu estimator with input μ: ",
                Ne,
            )
            time_to_coalescence = sum(
                (4 * Ne) / (i * (i - 1)) for i in range(2, total_seqs + 1)
            )
            self.coalescence_time = time_to_coalescence
            print(
                "Coalescence time using Effective Population size from input μ: ",
                time_to_coalescence,
            )
            self.mu = inp_mu

        return Ne

    def estimate_eff_pop_size_watterson(self, sequences):
        """
        Estimates the effective population size using the Watterson estimator

        Input: array of sequences
        Output: effective population size
        """
        total_seqs = len(sequences)
        # Mu - most often has the rate of 10e-4.  So we will use that value here
        mu = 0.00001
        # number of segregating sites
        K = 0
        # harmonic number n-1
        alpha_n = 0

        for i in range(1, len(sequences)):
            alpha_n += 1 / i

        threshold = 0.7
        tracking_dict = {}

        for i in range(len(sequences[0])):
            for seq in sequences:
                if seq[i] in tracking_dict:
                    tracking_dict[seq[i]] += 1
                else:
                    tracking_dict[seq[i]] = 1
            # determine if it is segregating site
            highest_char_cnt = -1
            for key in tracking_dict:
                if tracking_dict[key] > highest_char_cnt:
                    highest_char_cnt = tracking_dict[key]
            if highest_char_cnt / len(sequences) < threshold:
                K += 1
            tracking_dict = {}

        theta_w = K / alpha_n
        print("K: ", K)
        print("alpha_n: ", alpha_n)
        print("theta_w: ", theta_w)
        print("Default μ: ", mu)
        self.mu = mu
        Ne = theta_w / (4 * mu)
        print(
            "Effective Population size using Watterson estimator with default μ: ", Ne
        )

        time_to_coalescence = sum(
            (4 * Ne) / (i * (i - 1)) for i in range(2, total_seqs + 1)
        )
        self.coalescence_time = time_to_coalescence
        print(
            "Coalescence time using Effective Population size from default μ: ",
            time_to_coalescence,
        )

        inp_mu = self.promt_mutation_rate()

        if inp_mu is not None:
            Ne = theta_w / (4 * inp_mu)
            print(
                "Effective Population size using Watterson estimator with input μ: ", Ne
            )
            time_to_coalescence = sum(
                (4 * Ne) / (i * (i - 1)) for i in range(2, total_seqs + 1)
            )
            self.coalescence_time = time_to_coalescence
            print(
                "Coalescence time using Effective Population size from input μ: ",
                time_to_coalescence,
            )
            self.mu = inp_mu

        # apply μ correction
        # correction_coefficient = (K * len(sequences[0])) / time_to_coalescence
        # print("Correction Coefficient ω: ", correction_coefficient)

        # Ne = theta_w / (4*(mu*correction_coefficient))
        # print('Effective Population size using Watterson estimator with corrected μ: ', Ne)
        return Ne

    def estimate_eff_pop_size_watterson_no_input_mu(self, sequences):
        """
        Estimates the effective population size using the Watterson estimator

        Input: array of sequences
        Output: effective population size
        """
        total_seqs = len(sequences)
        # Mu - most often has the rate of 10e-4.  So we will use that value here
        # number of segregating sites
        K = 0
        # harmonic number n-1
        alpha_n = 0

        for i in range(1, len(sequences)):
            alpha_n += 1 / i

        threshold = 0.7
        tracking_dict = {}

        for i in range(len(sequences[0])):
            for seq in sequences:
                if seq[i] in tracking_dict:
                    tracking_dict[seq[i]] += 1
                else:
                    tracking_dict[seq[i]] = 1
            # determine if it is segregating site
            highest_char_cnt = -1
            for key in tracking_dict:
                if tracking_dict[key] > highest_char_cnt:
                    highest_char_cnt = tracking_dict[key]
            if highest_char_cnt / len(sequences) < threshold:
                K += 1
            tracking_dict = {}

        theta_w = K / alpha_n
        print("K: ", K)
        print("alpha_n: ", alpha_n)
        print("theta_w: ", theta_w)
        print("Default μ: ", self.mu)
        Ne = float(theta_w / (4 * self.mu))
        print(
            "Effective Population size using Watterson estimator with default μ: ", Ne
        )

        time_to_coalescence = sum(
            (4 * Ne) / (i * (i - 1)) for i in range(2, total_seqs + 1)
        )
        self.coalescence_time = time_to_coalescence
        print(
            "Coalescence time using Effective Population size from default μ: ",
            time_to_coalescence,
        )

        return Ne
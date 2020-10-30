from Bio.Align.Applications import ClustalOmegaCommandline, ClustalwCommandline
from Bio import SeqIO
from Bio import Phylo
from Bio import AlignIO
from Bio.Seq import Seq
from datetime import datetime
from util.file_util import FileStream
import numpy as np
import os
import random

file_tool = FileStream()

class SequenceUtil(object):
    """
    All functions related to working directly with sequences other than the motif
    finding functions which are contained in their own files
    """

    def align_sequences_w2_fasta(self, gene_name):
        '''
        Aligns sequences using the ClustalW2 executable

        input: gene_name --> to specify which fasta file to read from
        output:  Displays a phylogenetic tree
        '''

        print("Aligning sequences using ClustalW2...")
        time = self.get_time()
        in_file = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'fastas' + os.path.sep + '{0}.fasta'.format(gene_name)
        out_file = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'align' + os.path.sep + '{0}_out_{1}.aln'.format(gene_name, time)

        if os.name == 'nt':
            executable = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'executables' + os.path.sep + 'clustalw2.exe'
        else:
            executable = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'executables' + os.path.sep + 'clustalw2'
        
        clustalw_cline = ClustalwCommandline(executable, infile=in_file, outfile=out_file)
        stdout, stderr = clustalw_cline()
        align = AlignIO.read(out_file, "clustal")

        print("Finished aligning sequences.")
        print(align)

        self.display_phylo_tree(gene_name, in_file)

    
    def align_sequences_w2_file(self, filename):
        print("Aligning sequences using ClustalW2...")
        print("file name:", filename)
        time = self.get_time()
        in_file = os.getcwd() + os.path.sep + 'DNA' + os.path.sep + '{}.fasta'.format(filename)
        out_file = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'align' + os.path.sep + '{0}_out_{1}.aln'.format(filename, time)

        if os.name == 'nt':
            executable = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'executables' + os.path.sep + 'clustalw2.exe'
        else:
            executable = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'executables' + os.path.sep + 'clustalw2'

        clustalw_cline = ClustalwCommandline(executable, infile=in_file, outfile=out_file)
        stdout, stderr = clustalw_cline()
        align = AlignIO.read(out_file, "clustal")

        print("Finished aligning sequences.")
        print(align)

        self.display_phylo_tree(filename)

    
    def align_results_w2_single(self, generation_dict, generations, index, seq_id):
        print("Aligning results...")
        sequences_gens = [generation_dict[i][index] for i in range(generations)]
        file_tool.write_to_fasta_results(seq_id, sequences_gens)
        time = self.get_time()
        in_file = file_tool.most_recent_fasta_results()
        out_file = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'align' + os.path.sep + '{0}_{1}_aln'.format(seq_id.replace("|", ""), time)

        if os.name == 'nt':
            executable = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'executables' + os.path.sep + 'clustalw2.exe'
        else:
            executable = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'executables' + os.path.sep + 'clustalw2'

        clustalw_cline = ClustalwCommandline(executable, infile=in_file, outfile=out_file)
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
        out_file = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'align' + os.path.sep + 'post_sim_alignment_{}.aln'.format(time)

        if os.name == 'nt':
            executable = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'executables' + os.path.sep + 'clustalw2.exe'
        else:
            executable = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'executables' + os.path.sep + 'clustalw2'

        clustalw_cline = ClustalwCommandline(executable, infile=in_file, outfile=out_file)
        stdout, stderr = clustalw_cline()
        align = AlignIO.read(out_file, "clustal")

        print("Finished aligning sequences.")
        print(align)

        self.display_phylo_tree_results()


    def display_phylo_tree(self, filename):
        '''
        Reads from the dnd file produced by the ClustalW2 alignment
        and displays the phylogenetic tree based on the results of the alignment

        input: gene_name --> to specify the alignment file
               filename  --> specifies the path to the file
        output: Diplays phylogenetic tree using the Phylo biopython module
        '''

        print('Displaying Phylogenetic tree of sequences...')
        dnd_file = os.getcwd() + os.path.sep + 'DNA' + os.path.sep + '{}.dnd'.format(filename)
        tree = Phylo.read(dnd_file, 'newick')
        Phylo.draw(tree)

    
    def display_phylo_tree_results(self):
        print('Displaying Phylogenetic tree of sequences...')
        dnd_file = file_tool.most_recent_dnd_results()
        tree = Phylo.read(dnd_file, 'newick')
        Phylo.draw(tree)


    def calculate_conserved_regions(self):
        '''
        Naively determines conserved regions based on the alignment of the sequences

        input: no input but takes the most recent alignment file
        output: dictionary of each conserved region with the value consisting of a tuple with the start and end index, as well as the conserved DNA string
        '''
        #fetch aligned sequences from most recent alignment
        aligned_seqs = file_tool.read_from_alignment()
           
        print("Classifying conserved regions in the alignment")

        #declare tracking variables
        cur = 0
        start = 0
        chunk = ''
        chunks = {}
        chunk_count = 1

        #loop through each char in the sequences to determine whether they are conserved or not
        while cur < len(aligned_seqs[0]):
            column = [row[cur] for row in aligned_seqs]
            if self.all_same(column):
                chunk += aligned_seqs[0][cur]
            else:
                if len(chunk) > 1:
                    chunks[chunk_count] = (start, cur, chunk)
                    chunk_count += 1
                chunk = ''
                start = cur+1
            cur += 1
        print("Done classifying conserved regions")
        return chunks

    
    def check_for_cr(self, index, c_regions):
        '''
        checks if the nucleotide is a conserved region.  If it is, does not mutate
        '''
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

    
    def roll_duplication(self, time_elapsed):
       # alpha = .0000000000000001
       # prb = alpha*time_elapsed
       # roll = random.random()
       # if roll < prb:
        #    return True
        return False

    
    def roll_extinction(self, time_elapsed):
        pass

    
    def get_gc_content(self, seq):
        if len(seq) == 0:
            return 0
        return float((seq.count("G") + seq.count("C")) / len(seq))

    
    def get_nuc_count(self, seq):
        count_dict = {"A":0, "G":0, "T":0, "C": 0}
        for key in count_dict:
            count_dict[key] = seq.count(key)
        return count_dict


    def estimate_substitutions_generations(self, generation_dict, generations):
        first_gen = generation_dict[0]
        final_gen = generation_dict[generations - 1]

        #calculate divergence
        for i in range(len(first_gen)):
            score = 0
            origin_seq = first_gen[i]
            final_seq = final_gen[i]
            for j in range(len(origin_seq)):
                if origin_seq[j] == final_seq[j]:
                    score += 1
            print("Sequence {0} in the first and last generation are {1} percent similar".format(i, float(score/len(origin_seq)) * 100))


    def calculate_divergence_jc(self, generation_dict, generations):
        '''
        Calculates the K value (divergence) of the sequences post evolution for the 
        Jukes Cantor model
        '''

        first_gen = generation_dict[0]
        final_gen = generation_dict[generations - 1]

        for i in range(len(first_gen)):
            differences = 0
            origin_seq = first_gen[i]
            final_seq = final_gen[i]
            if len(origin_seq) == len(final_seq):
                for j in range(len(origin_seq)):
                    if origin_seq[j] != final_seq[j]:
                        differences += 1
            p = float(differences / len(origin_seq))
            k = -.75 * np.log(1 - 1.25*p)
            print("The K value for sequence {0} is {1}".format(i, k))

        
    def calculate_divergence_k2p(self, generation_dict, generations):
        '''
        Calculates the K value (divergence) of the sequences post evolution for the 
        Kimura 2P model
        '''

        transversions = {
            "A": ["C", "T"],
            "G": ["C", "T"],
            "C": ["A", "G"],
            "T": ["A", "G"],
            "-": []
        }

        transitions = {
            "A": ["G"],
            "G": ["A"],
            "C": ["T"],
            "T": ["C"],
            "-": []
        }

        first_gen = generation_dict[0]
        final_gen = generation_dict[generations - 1]

        for i in range(len(first_gen)):
            p = 0
            q = 0
            origin_seq = first_gen[i]
            final_seq = final_gen[i]
            if len(origin_seq) == len(final_seq):
                for j in range(len(origin_seq)):
                    char_origin = origin_seq[j]
                    char_final = final_seq[j]
                    if char_final in transversions[char_origin]:
                        q += 1
                    elif char_final in transitions[char_origin]:
                        p += 1
            tp = float(p / len(origin_seq))
            tv = float(q / len(origin_seq))
            k = -.5 * np.log(1-2*tp-tv) - .25 * np.log(1-2*tv)
            
            print("The K value for sequence {0} is {1}".format(i, k))


    def coalesce(self, sequences):
        '''
        Infers the ancestral sequence of a group of sequences generated from BLAST
        input:  array of DNA sequences
        output:  inferred ancestral sequence
        '''
        def merge(seq1, seq2):
            ancestor = ''
            for i in range(len(seq1)):
                s1 = seq1[i]
                s2 = seq2[i]

                if s1 == s2 and s1 != '-':
                    ancestor += s1
                elif s1 == '-' and s2 == '-':
                    ancestor += random.choice(['A', 'G', 'T', 'C'])
                elif s1 == '-' or s2 =='-':
                    ancestor += s1 if s1 != '-' else s2
                else:
                    nuc = random.choice([s1, s2])
                    ancestor += nuc
            return ancestor

        while len(sequences) > 1:

            rand1 = random.randint(0, len(sequences)-1)
            rand2 = random.randint(0, len(sequences)-1)

            while rand1 == rand2:
                rand1 = random.randint(0, len(sequences)-1)
                rand2 = random.randint(0, len(sequences)-1)
            
            seq1 = sequences[rand1]
            seq2 = sequences[rand2]

            temp1 = rand1 if rand1 < rand2 else rand2
            temp2 = rand2 if rand2 > rand1 else rand1

            sequences = sequences[:temp1] + sequences[temp1+1:temp2] + sequences[temp2+1:]

            new_seq = merge(seq1, seq2)

            sequences.append(new_seq)
        print("Ancestral sequence inferred: ", sequences[0])
        return sequences[0]




            



        
        

                    




                



from Bio.Align.Applications import ClustalOmegaCommandline, ClustalwCommandline
from Bio import SeqIO
from Bio import Phylo
from Bio import AlignIO
from datetime import datetime
from util.file_util import FileStream
import os

file_tool = FileStream()

class SequenceUtil(object):
    """
    All functions related to working directly with sequences other than the motif
    finding functions which are contained in their own files
    """

    def align_sequences_w2(self, gene_name):
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
            executable = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'executables' + os.path.sep + 'clustal-omega-1.2.3-macosx'
        
        clustalw_cline = ClustalwCommandline(executable, infile=in_file, outfile=out_file)
        stdout, stderr = clustalw_cline()
        align = AlignIO.read(out_file, "clustal")

        print("Finished aligning sequences.")
        print(align)

        self.display_phylo_tree(gene_name, in_file)


    def display_phylo_tree(self, gene_name, filename):
        '''
        Reads from the dnd file produced by the ClustalW2 alignment
        and displays the phylogenetic tree based on the results of the alignment

        input: gene_name --> to specify the alignment file
               filename  --> specifies the path to the file
        output: Diplays phylogenetic tree using the Phylo biopython module
        '''

        print('Displaying Phylogenetic tree of sequences...')
        dnd_file = filename.replace('{}.fasta'.format(gene_name), '{}.dnd'.format(gene_name))
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
        return chunks

    
    def all_same(self, items):
        return all(x == items[0] for x in items)
        

    def get_time(self):
        return datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
from Bio.Align.Applications import ClustalOmegaCommandline, ClustalwCommandline
from Bio import SeqIO
from Bio import Phylo
from Bio import AlignIO
import os
from datetime import datetime

class SequenceUtil(object):


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

    
    def get_time(self):
        return datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
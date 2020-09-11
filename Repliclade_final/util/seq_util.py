from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from Bio import Phylo
from Bio import AlignIO
import os
from datetime import datetime

class SequenceUtil(object):


    def align_sequences_omega(self, gene_name):
        print("Aligning sequences using ClustalOmega...")
        in_file = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'fastas' + os.path.sep + '{0}.fasta'.format(gene_name)
        out_file = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'align' + os.path.sep + '{0}_out_{1}.fasta'.format(gene_name, self.get_time())
        if os.name == 'nt':
            executable = ''
        else:
            executable = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'executables' + os.path.sep + 'clustal-omega-1.2.3-macosx'
        clustalomega_cline = ClustalOmegaCommandline(executable, outfmt = 'align', infile=in_file, outfile=out_file, verbose=True, auto=True)
        stdout, stderr = clustalomega_cline()
        align = AlignIO.read(out_file, "clustal")
        print("Finished aligning sequences.")
        self.display_phylo_tree(gene_name, out_file)


    def write_to_fasta(self, sequences, gene_name):
        DNA_dir = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'fastas' + os.path.sep + gene_name
        with open(DNA_dir+'.fasta', 'w') as open_file:
            SeqIO.write(sequences, open_file, "fasta")


    def display_phylo_tree(self, gene_name, filename):
        print('Displaying Phylogenetic tree of sequences...')
        alignment_path = filename
        tree = Phylo.read(alignment_path, 'newick')
        Phylo.draw()

    
    def get_time(self):
        return datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
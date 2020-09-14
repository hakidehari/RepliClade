from Bio import SeqIO
import os
import glob

class FileStream(object):

    def __init__(self, filename=None, filepath=None):
        self.filename = filename
        self.filepath = filepath
        self.cur_open_file = None

    
    def open_file(self):
        with open(self.filename) as f:
            self.cur_open_file = f

    def write_to_file(self, msg):
        if not self.cur_open_file:
            print("File not open for editing or not specified properly")
            return


    def most_recent_file(self):
        path_to_file = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'align' + os.path.sep + '*.aln'
        list_of_files = glob.glob(path_to_file)
        latest_file = max(list_of_files, key=os.path.getctime)
        print('Taking most recent file: {}'.format(latest_file))
        return latest_file


    def read_from_alignment(self):
        alignment_file = self.most_recent_file()
        aligned_seqs = []
        
        with open(alignment_file, "rU") as handle:
            for record in SeqIO.parse(handle, 'clustal'):
                aligned_seqs.append(str(record.seq))
        
        return aligned_seqs

    
    def read_from_fasta(self, gene_name):
        fasta_file = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'fastas' + os.path.sep + '{}.fasta'.format(gene_name)
        fasta_seqs = []

        with open(fasta_file, "rU") as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                fasta_seqs.append(str(record.seq))
                
        return fasta_seqs

    
    def write_to_fasta(self, sequences, gene_name):
        DNA_dir = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'fastas' + os.path.sep + gene_name
        with open(DNA_dir+'.fasta', 'w') as open_file:
            SeqIO.write(sequences, open_file, "fasta")
        
        
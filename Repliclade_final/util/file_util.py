from Bio import SeqIO
import os
import json
import glob
from Bio.Blast import NCBIXML
from datetime import datetime

class FileStream(object):


    def most_recent_file(self):
        '''
        Gets most recent file from the last alignment performed
        '''
        path_to_file = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'align' + os.path.sep + '*.aln'
        list_of_files = glob.glob(path_to_file)
        latest_file = max(list_of_files, key=os.path.getctime)
        print('Taking most recent file: {}'.format(latest_file))
        return latest_file


    def most_recent_file_results(self):
        '''
        Gets the most recent results json file
        '''
        path_to_file = os.getcwd() + os.path.sep + 'sim' + os.path.sep + 'results' + os.path.sep + '*.json'
        list_of_files = glob.glob(path_to_file)
        latest_file = max(list_of_files, key=os.path.getctime)
        return latest_file

    
    def most_recent_fasta_results(self):
        '''
        Gets the most recent results fasta file
        '''
        path_to_file = os.getcwd() + os.path.sep + 'sim' + os.path.sep + 'results' + os.path.sep + '*.fasta'
        list_of_files = glob.glob(path_to_file)
        latest_file = max(list_of_files, key=os.path.getctime)
        return latest_file

    
    def most_recent_dnd_results(self):
        '''
        Gets the most recent results dnd file
        '''
        path_to_file = os.getcwd() + os.path.sep + 'sim' + os.path.sep + 'results' + os.path.sep + '*.dnd'
        list_of_files = glob.glob(path_to_file)
        latest_file = max(list_of_files, key=os.path.getctime)
        return latest_file


    def read_from_alignment(self):
        '''
        Reads sequences from the most recent alignment file

        Outputs an array of aligned sequences
        '''
        alignment_file = self.most_recent_file()
        aligned_seqs = []
        
        with open(alignment_file, "rU") as handle:
            for record in SeqIO.parse(handle, 'clustal'):
                aligned_seqs.append(str(record.seq))
        
        return aligned_seqs

    
    def read_from_alignment_results(self):
        '''
        Reads sequences from the most recent alignment file

        Outputs an array of aligned sequences
        '''
        alignment_file = self.most_recent_file()
        aligned_seqs = []
        
        with open(alignment_file, "rU") as handle:
            for record in SeqIO.parse(handle, 'clustal'):
                aligned_seqs.append((str(record.seq), str(record.id)))
        
        return aligned_seqs

    
    def read_from_fasta(self, gene_name):
        '''
        reads from input fasta file before running the sequence against the BLAST algorithm
        '''
        fasta_file = os.getcwd() + os.path.sep + 'alignment' + os.path.sep + 'fastas' + os.path.sep + '{}.fasta'.format(gene_name)
        fasta_seqs = []

        with open(fasta_file, "rU") as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                fasta_seqs.append(str(record.seq))
                
        return fasta_seqs

    
    def read_from_blast_fasta(self, filename):
        '''
        Reads sequences from the fasta file generated by the blast function and used in the alignment
        '''
        fasta_file = os.getcwd() + os.path.sep + 'DNA' + os.path.sep + '{}.fasta'.format(filename)
        fasta_seqs = []

        with open(fasta_file, "rU") as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                fasta_seqs.append((str(record.seq), record.id))
                
        return fasta_seqs

    
    def read_from_blast(self, filename):
        '''
        Reads results from BLAST algorithm and returns them in an array of tuples
        '''
        blast_file_path = os.getcwd() + os.path.sep + 'DNA' + os.path.sep + '{}.xml'.format(filename)
        blast_seqs = []
        result_handle = open(blast_file_path)
        blast_records = NCBIXML.read(result_handle)
        for alignment in blast_records.alignments:
            if len(alignment.hsps) > 1:
                best_hsp = self.find_best_hsp(alignment.hsps)
                print("****Alignment****")
                print("sequence:", alignment.title)
                print("length:", alignment.length)
                print("e value:", best_hsp.expect)
                print(best_hsp.query[0:75] + "...")
                print(best_hsp.match[0:75] + "...")
                print(best_hsp.sbjct[0:75] + "...")
                blast_seqs.append((alignment.title, best_hsp.sbjct))
            else:
                for hsp in alignment.hsps:
                    print("****Alignment****")
                    print("sequence:", alignment.title)
                    print("length:", alignment.length)
                    print("e value:", hsp.expect)
                    print(hsp.query[0:75] + "...")
                    print(hsp.match[0:75] + "...")
                    print(hsp.sbjct[0:75] + "...")
                    blast_seqs.append((alignment.title, hsp.sbjct))
        return blast_seqs


    def find_best_hsp(self, hsps):
        '''
        Finds best matching hsp's for hsp's in a multiple hsp hit from BLAST
        '''
        max_score = -9999999
        best_hsp = None
        for hsp in hsps:
            if hsp.score > max_score:
                max_score = hsp.score
                best_hsp = hsp
        return best_hsp


    def write_to_fasta_sim_results(self, sequences, nodes, filename_original):
        from util.seq_util import SequenceUtil
        seq_obj = SequenceUtil()
        filename_original = filename_original.replace(".fasta", "").replace(".FASTA", "")
        filename = 'simulation_results_{}'.format(seq_obj.get_time())
        DNA_dir = os.getcwd() + os.path.sep + 'DNA' + os.path.sep
        with open(DNA_dir+'{}.fasta'.format(filename), 'w') as open_file:
            for seq in sequences:
                for node in nodes:
                    if node.get_sequence() == seq:
                        open_file.write('>{0}\n{1}\n'.format(filename_original + "_sim_" + str(node.num), seq))
        return filename

    
    def write_to_fasta_blast(self, sequences, filename):
        '''
        writes blast sequences to fasta
        '''
        DNA_dir = os.getcwd() + os.path.sep + 'DNA' + os.path.sep
        with open(DNA_dir+'{}.fasta'.format(filename), 'w') as open_file:
            for seq in sequences:
                open_file.write('>{0}\n{1}\n'.format(seq[0], seq[1]))

    
    def write_to_fasta_results(self, seq_id, sequences):
        dest_dir = os.getcwd() + os.path.sep + 'sim' + os.path.sep + 'results' + os.path.sep
        with open(dest_dir+'{}_results.fasta'.format(seq_id.replace("|", "")), "w") as open_file:
            for i in range(len(sequences)):
                open_file.write('>{0}Gen{1}\n{2}\n'.format(seq_id, i, sequences[i]))

    
    def write_to_fasta_results_multiple(self, seq_ids, sequences):
        dest_dir = os.getcwd() + os.path.sep + 'sim' + os.path.sep + 'results' + os.path.sep
        with open(dest_dir+'post_sim_results.fasta', 'w') as open_file:
            for i in range(len(sequences)):
                open_file.write('>{0}\n{1}\n'.format(seq_ids[i], sequences[i]))

    
    def log_simulation_to_json(self, generation_dict):
        '''
        Save results from each generation into a json file
        '''
        date_time = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        path_to_file = os.getcwd() + os.path.sep + 'sim' + os.path.sep + 'results' + os.path.sep
        filename = 'sim_results_{}.json'.format(date_time)
        with open(path_to_file + filename, 'w') as open_file:
            open_file.write(json.dumps(generation_dict))
        print("Results written to {}".format(path_to_file + filename))

    
    def log_simulation_output_to_json(self, output_dict):
        date_time = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        path_to_file = os.getcwd() + os.path.sep + 'sim' + os.path.sep + 'results' + os.path.sep
        filename = 'sim_output_{}.json'.format(date_time)
        with open(path_to_file + filename, 'w') as open_file:
            open_file.write(json.dumps(output_dict))
        print("Results written to {}".format(path_to_file + filename))


    def read_simulation_results(self):
        latest_file_path = self.most_recent_file_results()
        generation_dict = json.loads(latest_file_path)
        return generation_dict


        
        
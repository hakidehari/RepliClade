import sys
import os
from config import config
from DNA import *
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import base64


class GenBankConnector(object):
    def __init__(self):
        """initializes Genbank object to connect to Genbank"""
        Entrez.api_key = base64.b64decode(config.genbank_auth["api_key"]).decode(
            "utf-8"
        )
        Entrez.email = base64.b64decode(config.genbank_auth["email"]).decode("utf-8")

    def fetch_sequences(self, sequences):
        """
        Fetches sequences from GenBank given an array of GenBank Id's or singular Id

        input:  array or string
        returns: array of sequence objects
        """

        print("Fetching sequence from GenBank...")
        seq_array = []
        if isinstance(sequences, list):
            for seq_id in sequences:
                handle = Entrez.efetch(
                    db="nucleotide", id=seq_id, rettype="gb", retmode="text"
                )
                record = SeqIO.read(handle, "genbank")
                seq_array.append(record)
        else:
            handle = Entrez.efetch(
                db="nucleotide", id=sequences, rettype="gb", retmode="text"
            )
            record = SeqIO.read(handle, "genbank")
            seq_array.append(record)
        print("Finished fetching sequences from GenBank.")
        return seq_array

    def run_ncbi_blast(self, gene):
        print("Executing BLAST on the sequences")
        fasta_path = (
            os.getcwd()
            + os.path.sep
            + "alignment"
            + os.path.sep
            + "fastas"
            + os.path.sep
        )
        record = SeqIO.read(open(fasta_path + "{}.fasta".format(gene)), format="fasta")
        result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))
        if os.name == "nt":
            save_file_path = os.getcwd() + "\\alignment\\blast\\{}.xml".format(gene)
        else:
            save_file_path = os.getcwd() + "/alignment/blast/{}.xml".format(gene)
        save_file = open(save_file_path, "w")
        save_file.write(result_handle.read())
        save_file.close()
        result_handle.close()

    def run_ncbi_blast_input_file(self, filename):
        print("Executing BLAST on the sequences")
        fasta_path = os.getcwd() + os.path.sep + "DNA" + os.path.sep + filename
        record = SeqIO.read(open(fasta_path), format="fasta")
        result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))
        if os.name == "nt":
            save_file_path = os.getcwd() + "\\DNA\\{}.xml".format(filename)
        else:
            save_file_path = os.getcwd() + "/DNA/{}.xml".format(filename)
        save_file = open(save_file_path, "w")
        save_file.write(result_handle.read())
        save_file.close()
        result_handle.close()

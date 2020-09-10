import sys
from config import config
from DNA import *
from Bio import Entrez
from Bio import SeqIO
import base64


class GenBankConnector(object):

    def __init__(self):
        '''initializes Genbank object to connect to Genbank'''
        Entrez.api_key = base64.b64decode(config.genbank_auth['api_key']).decode('utf-8')
        Entrez.email = base64.b64decode(config.genbank_auth['email']).decode('utf-8')


    def fetch_sequences(self, sequences):
        '''
        Fetches sequences from GenBank given an array of GenBank Id's or singular Id

        input:  array or string
        returns: array of sequence objects
        '''

        print("Fetching sequence from GenBank...")
        seq_array = []
        if isinstance(sequences, list):
            for seq_id in sequences:
                handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
                record = SeqIO.read(handle, "genbank")
                seq_array.append(record.seq)
        else:
            handle = Entrez.efetch(db="nucleotide", id=sequencea, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            seq_array.append(record.seq)
        print("Finished fetching sequences from GenBank.")
        return seq_array


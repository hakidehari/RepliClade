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

    def fetch_sequences(self, filename):
        pass


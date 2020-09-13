

class DetMotifFinding:

    def __init__(self, size=8, seqs=None):
        self.motif_size = size
        if seqs is not None:
            self.seqs = seqs
            self.alphabet = seqs[0].alphabet()
        else:
            self.seqs = []

    
    def __len__(self):
        return len(self.seqs)


    def __getitem__(self, n):
        return self.seqs[n]

    
    def seq_size(self, i):
        return len(self.seqs[i])
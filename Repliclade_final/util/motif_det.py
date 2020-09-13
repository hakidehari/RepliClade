

class DetMotifFinding:

    def __init__(self, size=8, seqs=None):
        self.motif_size = size
        if seqs is not None:
            self.seqs = seqs
            self.alphabet = self.get_alphabet(seqs[0])
        else:
            self.seqs = []


    def get_alphabet(self, s):
        res = set()
        for char in s:
            res.add(char)
        return list(res)


    def __len__(self):
        return len(self.seqs)


    def __getitem__(self, n):
        return self.seqs[n]

    
    def seq_size(self, i):
        return len(self.seqs[i])


    def create_motifs_from_indexes(self, indexes):
        pseqs = []
        res = [[0] * self.motif_size for char in range(len(self.alphabet))]

        for i, ind in enumerate(indexes):
            subseq = self.seqs[i][ind:(ind + self.motif_size)]
            for i in range(self.motif_size):
                for k in range(len(self.alphabet)):
                    if subseq[i] == self.alphabet[k]:
                        res[k][i] = res[k][i] + 1
        
        return res
    

    def score(self, s):
        score = 0
        mat = self.create_motifs_from_indexes(s)
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for i in range(1, len(mat)):
                if mat[i][j] > maxcol:
                    maxcol = mat[i][j]
            score += maxcol
        
        return score

    
    def score_multiplicative(self, s):
        score = 1.0
        mat = self.create_motifs_from_indexes(s)
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for i in range(1, len(mat)):
                if mat[i][j] > maxcol:
                    maxcol = mat[i][j]
            score *= maxcol
        return score

    
    def next_solution(self, s):
        next_sol = [0]*len(s)
        pos = len(s) - 1
        while pos >= 0 and s[pos] == self.seq_size(pos) - self.motif_size:
            pos -= 1
        if pos < 0:
            next_sol = None
        else:
            for i in range(pos):
                next_sol[i] = s[i]
            next_sol[pos] = s[pos] + 1
            for i in range(pos+1, len(s)):
                next_sol[i] = 0
        return next_sol


    def exhaustive_search(self):
        best_score = -1
        res = []
        s = [0]*len(self.seqs)
        while s is not None:
            sc = self.score(s)
            if sc > best_score:
                best_score = sc
                res = s
            s = self.next_solution(s)
        return res

    
    def next_vertex(self, s):
        res = []
        if len(s) < len(self.seqs):
            for i in range(len(s)):
                res.append(s[i])
            res.append(0)
        else:
            pos = len(s) - 1
            while pos >= 0 and s[pos] == self.seq_size(pos) - self.motif_size:
                pos -= 1
            if pos < 0:
                res = None
            else:
                for i in range(pos):
                    res.append(s[i])
                res.append(s[pos] + 1)
        return res

    
    def bypass(self, s):
        res = []
        pos = len(s) - 1
        while pos >= 0 and s[pos] == self.seq_size(pos) - self.motif_size:
            pos -= 1
        if pos < 0:
            res = None
        else:
            for i in range(pos):
                res.append(s[i])
            res.append(s[pos] + 1)
        return res

    
    def branch_and_bound(self):
        best_score = -1
        best_motif = None
        size = len(self.seqs)
        s = [0]*size
        while s is not None:
            if len(s) < size:
                optimum_score = self.score(s) + (size-len(s)) * self.motif_size
                if optimum_score < best_score:
                    s = self.bypass(s)
                else:
                    s = self.next_vertex(s)
            else:
                sc = self.score(s)
                if sc > best_score:
                    best_score = sc
                    best_motif = s
                s = self.next_vertex(s)
        return best_motif

    
    def heuristic_consensus(self):
        res = [0] * len(self.seqs)
        max_score = -1
        partial = [0, 0]
        for i in range(self.seq_size(0) - self.motif_size):
            for j in range(self.seq_size(1) - self.motif_size):
                partial[0] = i
                partial[1] = j
                sc = self.score(partial)
                if sc > max_score:
                    max_score = sc
                    res[0] = i
                    res[1] = j
        for k in range(2, len(self.seqs)):
            partial = [0]*(k+1)
            for j in range(k):
                partial[j] = res[j]
            max_score = -1
            for i in range(self.seq_size(k) - self.motif_size):
                partial[k] = i
                sc = self.score(partial)
                if sc > max_score:
                    max_score = sc
                    res[k] = i
        return res
            

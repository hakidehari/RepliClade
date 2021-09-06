class MotifPrb:
    def __init__(self, seqs=[], pwm=[], alphabet=None):
        if seqs:
            self.size = len(seqs[0])
            self.seqs = seqs
            self.alphabet = ["A", "G", "T", "C", "-"]
            self.do_counts()
            self.create_pwm()
        else:
            self.pwm = pwm
            self.size = len(pwm[0])
            self.alphabet = alphabet

    def __len__(self):
        return self.size

    def do_counts(self):
        self.counts = self.create_matrix_zeros(len(self.alphabet), self.size)
        for s in self.seqs:
            for i in range(self.size):
                lin = self.alphabet.index(s[i])
                self.counts[lin][i] += 1

    def create_pwm(self):
        if self.counts is None:
            self.do_counts()
        self.pwm = self.create_matrix_zeros(len(self.alphabet), self.size)
        for i in range(len(self.alphabet)):
            for j in range(self.size):
                self.pwm[i][j] = float(self.counts[i][j] / len(self.seqs))

    def consensus(self):
        res = ""
        for j in range(self.size):
            maxcol = self.counts[0][j]
            maxcoli = 0
            for i in range(1, len(self.alphabet)):
                if self.counts[i][j] > maxcol:
                    maxcol = self.counts[i][j]
                    maxcoli = i
            res += self.alphabet[maxcoli]
        return res

    def masked_consenus(self):
        res = ""
        for j in range(self.size):
            maxcol = self.counts[0][j]
            maxcoli = 0
            for i in range(1, len(self.alphabet)):
                if self.counts[i][j] > maxcol:
                    maxcol = self.counts[i][j]
                    maxcoli = i
            if maxcol > len(self.seqs / 2):
                res += self.alphabet[maxcoli]
            else:
                res += "-"
        return res

    def probability_sequence(self, seq):
        res = 1.0
        for i in range(self.size):
            lin = self.alphabet.index(seq[i])
            res *= self.pwm[lin][i]
        return res

    def probability_all_positions(self, seq):
        res = []
        for k in range(len(seq) - self.size + 1):
            res.append(self.probability_sequence(seq))
        return res

    def most_probable_sequence(self, seq):
        maximum = -1.0
        maxind = -1
        for k in range(len(seq) - self.size):
            p = self.probability_sequence(seq[k : k + self.size])
            if p > maximum:
                maximum = p
                maxind = k
        return maxind

    def create_motif(self, seqs):
        l = []
        for s in seqs:
            ind = self.most_probable_sequence(s)
            subseq = s[ind : ind + self.size]
            l.append(subseq)
        return PrbMotifFinding(l)

    def get_alphabet(self, s):
        res = set()
        for seq in s:
            for char in seq:
                res.add(char)
        return list(res)

    def create_matrix_zeros(self, nrows, ncols):
        return [[0] * ncols for i in range(nrows)]

    def print_matrix(self, mat):
        for i in range(len(mat)):
            print(mat[i])

def print_tree(root):

    if root == None:
        return
   
    # Standard level order traversal code
    # using queue
    q = []  # Create a queue
    q.append(root); # Enqueue root 
    while len(q) != 0:
     
        n = len(q)
  
        # If this node has children
        while n > 0:
         
            # Dequeue an item from queue and print it
            p = q[0]
            q.pop(0)
            print(p.num, end=' ')
   
            # Enqueue all children of the dequeued item
            for i in range(len(p.get_children())):
             
                q.append(p.get_children()[i])
            n -= 1
   
        print()


class TreeNode(object):


    def __init__(self, seq, root=False, num=None):
        self.__sequence = seq
        self.__children = []
        self.__root = root
        self.num = num

    
    def get_sequence(self):
        return self.__sequence

    
    def get_children(self):
        return self.__children

    
    def add_child(self, node):
        self.__children.append(node)

    
    def set_sequence(self, seq):
        self.__sequence = seq

    
    def is_root(self):
        return self.__root

    
    def print_contents(self):
        print("Sequence {0}\n{1}".format(self.num, self.__sequence))


class Phylogenize(object):

    def __init__(self, seq=None):
        self.__sequences = seq

    
    def calculate_distance_jc(self):
        '''
        Calculates pairwise distances for the Jukes Cantor Model
        '''

        import numpy as np

        dist_matrix = [[0] * len(self.__sequences) for seq in self.__sequences]
        for i in range(len(self.__sequences)):
            for j in range(len(self.__sequences)):
                if i == j:
                    pass
                else:
                    differences = 0
                    seq1 = self.__sequences[i]
                    seq2 = self.__sequences[j]
                    for k in range(len(seq1)):
                        if seq1[k] != seq2[k]:
                            differences += 1
                    p = float(differences / len(seq1))
                    k = -.75 * np.log(1-1.25*p)
                    dist_matrix[i][j] = k

        return dist_matrix

    
    def calculate_distance_k2p(self):
        '''
        Calculates pairwise distances for the Kimura 2P Model
        '''

        import numpy as np

        transversions = {
            "A": ["C", "T"],
            "G": ["C", "T"],
            "C": ["A", "G"],
            "T": ["A", "G"],
            "-": []
        }

        transitions = {
            "A": ["G"],
            "G": ["A"],
            "C": ["T"],
            "T": ["C"],
            "-": []
        }

        dist_matrix = [[0] * len(self.__sequences) for seq in self.__sequences]
        for i in range(len(self.__sequences)):
            for j in range(len(self.__sequences)):
                if i == j:
                    pass
                else:
                    seq1 = self.__sequences[i]
                    seq2 = self.__sequences[j]
                    q = 0
                    p = 0
                    for k in range(len(seq1)):
                        char1 = seq1[k]
                        char2 = seq2[k]
                        if char2 in transversions[char1] or char1 in transversions[char2]:
                            q += 1
                        if char2 in transitions[char1] or char1 in transitions[char2]:
                            p += 1
                    tp = float(p / len(seq1))
                    tv = float(q / len(seq1))
                    k = -.5 * np.log(1-2*tp-tv) - .25 * np.log(1-2*tv)
                    dist_matrix[i][j] = k
        
        return dist_matrix

    
    def biopython_calc_distances_upgma_nj(self):
        from Bio.Phylo.TreeConstruction import DistanceCalculator
        from Bio import AlignIO
        from util.file_util import FileStream

        file_tool = FileStream()
        #read most recent alignment file
        alignment_file = file_tool.most_recent_file()
        aln = AlignIO.read(alignment_file, 'clustal')
        #calculator
        calculator = DistanceCalculator('identity')
        #calculate distance matrix
        dm = calculator.get_distance(aln)

        return calculator, dm

    
    def build_tree_upgma_nj(self, calculator, dm, method):
        from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
        from Bio import Phylo
        #neighbor joining = 'nj', UPGMA = 'upgma'
        TREE_LIMIT = 500
        trees = []

        
        constructor = DistanceTreeConstructor()
        if method == 'nj':
            tree = constructor.nj(dm)
        else:
            tree = constructor.upgma(dm)

        print(tree)
            

        Phylo.draw(tree)

        
    
    def prompt_tree_builder(self):
        choices = ['upgma', 'parsimony', 'nj', 'ml']
        print("How would you like to build the phylogenetic tree of the simulation?")
        method = input("You can select nj=neighbor joining, upgma=UPGMA, parsimony=Maximum Parsimony, ml=Maximum Likelihood: ")
        while method not in choices:
            method = input("Invalid choice, please choose one of the options above (nj, upgma, parsimony, or ml): ")
        return method

    
    def build_tree_parsimony(self):
        from Bio import AlignIO
        from Bio.Phylo.TreeConstruction import ParsimonyScorer, NNITreeSearcher, ParsimonyTreeConstructor
        from util.file_util import FileStream
        from Bio import Phylo
        


        file_tool = FileStream()
        #read from most recent alignment file
        alignment_file = file_tool.most_recent_file()
        aln = AlignIO.read(alignment_file, 'clustal')
        #instantiate parsimony scorer and NNI Tree Searcher

        trees = []
        
        scorer = ParsimonyScorer()
        searcher = NNITreeSearcher(scorer)
        constructor = ParsimonyTreeConstructor(searcher)
        #build parsimony tree
        pars_tree = constructor.build_tree(aln)
        print(pars_tree)

        Phylo.draw(pars_tree)

    
    def maximum_likelihood(self):
        ############IN PROGRESS###############
        import os
        from Bio import AlignIO
        from util.file_util import FileStream
        from Bio import Phylo
        from Bio.Phylo.Applications import PhymlCommandline

        print(os.getcwd())

        phyml_ex_path = os.getcwd() + os.sep + 'alignment' + os.sep + 'executables' + os.sep
        if os.name == 'nt':
            phyml_ex_path = phyml_ex_path + 'PhyML-3.1_win32.exe'
        else:
            phyml_ex_path = phyml_ex_path + 'PhyML-3.1_macOS-MountainLion'
        


        file_tool = FileStream()
        alignment_file = file_tool.most_recent_file()
        converted_file = os.getcwd()+os.sep+'alignment'+os.sep+'align'+os.sep+'converted_phyl.phy'
        aln = AlignIO.convert(alignment_file,
                             'clustal',
                             converted_file,
                             'phylip-relaxed')
    
        phyml = PhymlCommandline(cmd=phyml_ex_path, input=converted_file)
        out_log, err_log = phyml()
        dnd_file = converted_file + '_phyml_tree.txt'
        tree = Phylo.read(dnd_file, 'newick')
        Phylo.draw(tree)















        

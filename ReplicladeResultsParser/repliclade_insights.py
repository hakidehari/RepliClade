import os
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin


def get_all_file_names():
    '''
    Gets path to tree files and the filenames
    '''
    tree_file_path = os.path.dirname(os.getcwd()) + os.path.sep + 'Repliclade_final' + os.path.sep + 'sim' + os.path.sep + 'results' + os.path.sep + 'trees'
    files = []
    for (dirpath, dirnames, filenames) in os.walk(tree_file_path):
        files.extend(filenames)
        break
    return files, tree_file_path


def read_phylo_trees(filenames, path):
    '''
    Parses the trees from the files
    '''
    all_trees = []
    for file in filenames:
        trees = Phylo.parse(path + os.path.sep + file, 'phyloxml')
        arr = []
        for tree in trees:
            arr.append(tree)
        all_trees.append(arr)
    return all_trees


def clean_true_trees(trees):
    '''
    Cleans the names of the true trees to reflect names of other trees
    '''
    for treefile in trees:
        for clade in treefile[4].find_clades():
            clade.name = 'cytochrome_b_1_sim_{}'.format(clade.name)

    return trees
        

def calculate_tree_similarity(trees):
    #nj upgma parsimony ml
    scores = {}
    for file in trees:
        score = {}
        true_tree = file[4]
        for i in range(4):
            other_tree = file[i]
            true_terminals = [clade.name for clade in TreeMixin.get_terminals(true_tree)]
            other_terminals = [clade.name for clade in TreeMixin.get_terminals(other_tree)]
            true_nonterminals = [clade.name for clade in TreeMixin.get_nonterminals(true_tree)]
            other_nonterminals = [clade.name for clade in TreeMixin.get_nonterminals(other_tree)]

            count_true_terminals = TreeMixin.count_terminals(true_tree)
            count_other_terminals = TreeMixin.count_terminals(other_tree)
            return

def clean_helper(clade):
    pass
    




if __name__ == '__main__':
    filenames, path = get_all_file_names()
    trees = read_phylo_trees(filenames, path)

    cleaned_trees = clean_true_trees(trees)
    

    calculate_tree_similarity(trees)
    

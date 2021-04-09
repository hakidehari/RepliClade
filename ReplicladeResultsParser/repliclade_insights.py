import os
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin


def get_all_file_names(model, gene):
    '''
    Gets path to tree files and the filenames
    '''
    tree_file_path = os.path.dirname(os.getcwd()) + os.path.sep + 'Repliclade_final' + os.path.sep + 'sim' + os.path.sep + 'results' + os.path.sep + 'trees'
    files = []
    for (dirpath, dirnames, filenames) in os.walk(tree_file_path):
        files.extend(filenames)
        break
    return [file for file in files if model in file and gene in file], tree_file_path


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


def clean_true_trees(trees, gene):
    '''
    Cleans the names of the true trees to reflect names of other trees
    '''
    dic = {'cytochrome': 'cytochrome_b_1_sim_', 'TPF': 'TPF_sim_', 'AGP': 'AGP_sim_'}
    for treefile in trees:
        for clade in treefile[4].find_clades():
            clade.name = f'{dic[gene]}{clade.name}'
    return trees
        

def calculate_tree_similarity(trees):
    #nj upgma parsimony ml
    scores = {}
    index = 0
    for file in trees:
        score = {}
        true_tree = file[4]
        for i in range(4):
            other_tree = file[i]
            true_terminals = [clade.name for clade in TreeMixin.get_terminals(true_tree)]
            other_terminals = [clade.name for clade in TreeMixin.get_terminals(other_tree) if clade.name is not None and "Inner" not in clade.name]
            if TreeMixin.get_nonterminals(true_tree) is not None and TreeMixin.get_nonterminals(other_tree) is not None:
                true_nonterminals = [clade.name for clade in TreeMixin.get_nonterminals(true_tree)]
                other_nonterminals = [clade.name for clade in TreeMixin.get_nonterminals(other_tree) if clade.name is not None and "Inner" not in clade.name]
            else:
                true_nonterminals = []
                other_nonterminals = []

            terminal_nonterminal_score = compare_terminals(true_terminals, other_terminals, true_nonterminals, other_nonterminals, true_tree, other_tree)
            true_tree_clade_count = 0
            for clade in true_tree.find_clades():
                true_tree_clade_count += 1
            
            true_total_branch_len = TreeMixin.total_branch_length(true_tree)
            other_total_branch_len = TreeMixin.total_branch_length(other_tree)

            branch_len_score = abs(other_total_branch_len - true_total_branch_len)

            #parse
            score[i] = terminal_nonterminal_score# + branch_len_score
        
        scores[index] = score
        index += 1
    #print(scores)
    return scores


def compare_terminals(true_terminals, other_terminals, true_nonterminals, other_nonterminals, true_tree, other_tree):

    true_terminals = set(true_terminals)
    other_terminals = set(other_terminals)
    len_intersecting_terminals = len(true_terminals.intersection(other_terminals))
    total_terminals = len(true_terminals) + len(other_terminals) - len_intersecting_terminals
    terminals_final_score = len_intersecting_terminals / total_terminals
    
    true_nonterminals = set(true_nonterminals)
    other_nonterminals = set(other_nonterminals)
    len_intersecting_nonterminals = len(true_nonterminals.intersection(other_nonterminals))
    total_nonterminals = len(true_nonterminals) + len(other_nonterminals) - len_intersecting_nonterminals
    nonterminals_final_score = len_intersecting_terminals / total_nonterminals

    preterminal_count_true = 0
    preterminal_count_other = 0

    for clade in true_tree.find_clades():
        if TreeMixin.is_preterminal(clade):
            preterminal_count_true += 1
    for clade in other_tree.find_clades():
        if TreeMixin.is_preterminal(clade):
            preterminal_count_other += 1

    preterminal_score = 1 if preterminal_count_true == preterminal_count_other else abs(preterminal_count_other - preterminal_count_true)
    score = terminals_final_score + nonterminals_final_score / preterminal_score
    return score


def evaluate_scores(scores):

    tree_method_dict = {0:'nj', 1:'upgma', 2:'parsimony', 3:'ml'}

    total_scores = {}
    
    for sim_index in scores:
        for score in scores[sim_index]:
            if tree_method_dict[score] not in total_scores:
                total_scores[tree_method_dict[score]] = scores[sim_index][score]
            else:
                total_scores[tree_method_dict[score]] += scores[sim_index][score]
    #print(total_scores)
    return total_scores
    




if __name__ == '__main__':
    models = ['jukescantor', 'kimura', 'hasegawa', 'felsenstein']
    genes = ['cytochrome', 'TPF', 'AGP']
    all_scores = {}
    for gene in genes:
        print(gene)
        for model in models:
            filenames, path = get_all_file_names(model, gene)
            trees = read_phylo_trees(filenames, path)

            cleaned_trees = clean_true_trees(trees, gene)
            
            scores = calculate_tree_similarity(trees)

            all_scores[model] = evaluate_scores(scores)
        
        print(all_scores)
        all_scores = {}

    

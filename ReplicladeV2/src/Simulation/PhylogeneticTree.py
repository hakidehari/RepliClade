from Sequence import Sequence

class PhylogeneticTree(object):


    def __init__(self, root):
        self.root = root
        self.total_nodes = 0
        self.total_extinct_nodes = 0
        self.tree_depth = 0

    
    def display_tree(self):
        pass

    def get_max_tree_depth(self):
        pass


class TreeNode(object):


    def __init__(self, name=None, sequence=None, children=[]):
        self.name = name
        self.sequence = sequence
        self.children = children
        self.is_extinct = False
        self.gene_annotations = {}

    
    def add_children(self, child):
        '''Adds child to a node in the tree'''
        if not self.is_extinct:
            self.children.append(child)
        else:
            print("Cannot add a child to an extinct node.")

    
    def make_node_extinct(self):
        '''Makes a node extinct'''
        self.is_extinct = True


    def add_gene_annotation(self, start_index, end_index, gene_annotation, name=None):
        '''Adds a description/annotation to a specific range in the genome'''
        if self.is_extinct:
            print("Unable to add Gene Annotation.  Node is extinct")
            return

        for key in self.gene_annotations:
            annotation_tuple = self.gene_annotations[key]
            if annotation_tuple[0] in range(start_index, end_index + 1) or annotation_tuple[1] in range(start_index, end_index + 1):
                print("Unable to add Gene annotation.  Annotation already exists within the range specified")
                return

        if name is None:
            new_annotation = (start_index, end_index, gene_annotation)
        else:
            new_annotation = (start_index, end_index, gene_annotation, name)
            
        self.gene_annotations[len(self.gene_annotations)] = new_annotation


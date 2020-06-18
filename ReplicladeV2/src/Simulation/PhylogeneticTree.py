

class PhylogeneticTree(object):


    def __init__(self, name, root):
        self.root = root


class TreeNode(object):


    def __init__(self, name, children):
        self.name = name
        self.children = children
        self.is_extinct = False

    
    def add_child(self, child):
        self.children.append(child)

    
    def make_node_extinct(self):
        self.is_extinct = True
        
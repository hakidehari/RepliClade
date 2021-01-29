class TreeNode(object):

    
    def __init__(self, seq):
        self.seq = seq
        self.children = []
        self.is_extinct = False
        self.parent = None


    def add_children(self, sequences):
        if not self.is_extinct:
            if isinstance(sequences, list):
                self.children = [TreeNode(seq) for seq in sequences]
                for child in self.children
            else:
                print("No List provided for TreeNode.add_children function.")

    
    def make_extinct(self):
        self.is_extinct = True

    
    def get_parent(self):
        return self.parent
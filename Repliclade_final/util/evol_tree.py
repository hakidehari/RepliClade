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

    def __init__(self, seq):
        self.__sequences = seq




        

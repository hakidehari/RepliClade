
class Lineage(object):

    def __init__(self, seq):
        self.__sequence = seq
        self.__neighbors = []

    
    def add_neighbor(self, neighbor, relation, propogate):
        new_neighbor = Lineage(neighbor)
        if propogate:
            if relation == 'parent':
                new_neighbor.add_neighbor(self.__sequence, 'child', False)
            else:
                new_neighbor.add_neighbor(self.__sequence, 'parent', False)
        self.__neighbors.append(
            {
                'sequence': Lineage(neighbor), 
                'relation': relation
            }
        )

    
    def update_sequence(self, seq):
        self.__sequence = seq
        

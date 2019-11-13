import random
from Sequence import Sequence


class ProbabilityModels:

    def possibleIndelInsertion(self):
        '''Checks Whether an indel should happen or not'''
        r = random.random()
        if r < .00002:
            return True
        return False


    def influenzaMutate(self, sequence):
        '''Calculates the probability of a mutation per cell infection of a strain of the Influenze A virus'''
        for i in range (0, len(sequence.sequence)):
            if random.random() <= .000045:
                sequence = Sequence(sequence.modify(i))
        return sequence


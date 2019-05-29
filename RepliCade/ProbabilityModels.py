import random
from Sequence import Sequence


'''Checks Whether an indel should happen or not'''
def possibleIndelInsertion():
    r = random.random()
    if r < .00002:
        return True
    return False


'''Will run through a probability model and return whether the sequence should duplicate or not'''
def possibleDuplication():
    #for now returns true until a probability model is implemented
    r = random.random()
    if r < .00002:
        return True
    return False


'''Will run through a probability model and return whether a mutation will occur during duplication'''
def possibleMutationDuringDuplication(sequence):
    #for now assign a random probability until I can get a good probability model in here
    r = random.random()
    if r < 0.0002:
        sequence = Sequence(sequence.randomModify())
    return sequence


'''will run through a probability model and return whether the sequence should mutate or not'''
def possibleMutation():
    r = random.random()
    if r < .0005:
        return True
    return False

'''Calculates the probability of a mutation per cell infection of a strain of the Influenze A virus'''
def influenzaMutate(sequence):
    for i in range (0, len(sequence.sequence)):
        if random.random() <= .000045:
            print("hell yeah")
            sequence = Sequence(sequence.modify(i))
    return sequence


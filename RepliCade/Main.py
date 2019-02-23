"""

This is the Main file for Repliclade where the program will be run
All dependencies and classes will be compiled here and ran accordingly

"""


import random
from Generator import Generator
from Sequence import Sequence



'''returns the user specified DNA length and simulation time (in iterations)'''
def simulationParameters(DNALength,simulationTime):
    return [DNALength, simulationTime]


'''Will run through a probability model and return whether the sequence should duplicate or not'''
def possibleDuplication():
    #for now returns true until a probability model is implemented
    r = random.random()
    if r < .0002:
        return True
    return False


'''Will run through a probability model and return whether a mutation will occur during duplication'''
def possibleMutationDuringDuplication(sequence):
    #for now assign a random probability until I can get a good probability model in here
    r = random.random()
    if r < 0.5:
        sequence = Sequence(sequence.randomModify())
    return sequence


'''will run through a probability model and return whether the sequence should mutate or not'''
def possibleMutation():
    r = random.random()
    if r < .2:
        return True
    return False


'''
    Run simulation given two user specified parameters returned from the simulationParameters() function
    Currently runs smoothly with 10k iterations but slows down tremendously after 50k iterations
'''
def runSimulation():
    #vals[0] is the DNA length and vals[1] is the simulation time in (units here)
    vals = simulationParameters(20, 10000)
    seqLen = vals[0]
    runTime = vals[1]
    sequences = []
    a = Generator.generateSequenceGCContent(Generator, seqLen, 0.6)
    ancestor = Sequence(a)
    sequences.append(ancestor)

    #once the generation parameter runs out of steam, break from the while loop
    while runTime > 0:
        for sequence in sequences:
            #check if duplication should happen
            if possibleDuplication():
                #append with a possible error in duplication
                sequences.append(possibleMutationDuringDuplication(sequence))
            else:
                #else check the probability of a mutation happening
                if possibleMutation():
                    #if so, modify it and save it to the sequence list
                    sequence = Sequence(sequence.randomModify())
        runTime -= 1

    print(len(sequences))
    for sequence in sequences:
        print(sequence.sequence)
    return sequences



from Generator import Generator
from Sequence import Sequence



'''returns the user specified DNA length and simulation time (in generations)'''
def simulationParameters(DNALength,simulationTime):
    return [DNALength, simulationTime]

'''Will run through a probability model and return whether the sequence should duplicate or not'''
def possibleDuplication():
    #for now returns true until a probability model is implemented
    return True

'''Will run through a probability model and return whether a mutation will occur during duplication'''
def possibleMutationDuringDuplication(sequence):
    #for now true until I can get a good probability model in here
    sequence.randomModify()
    return sequence

'''will run through a probability model and return whether the sequence should mutate or not'''
def possibleMutation():
    return True


'''Run simulation given two user specified parameters returned from the simulationParameters() function'''
def runSimulation():
    #vals[0] is the DNA length and vals[1] is the simulation time in generations
    vals = simulationParameters(1000, 100)
    sequences = []
    ancestor = Sequence(Generator.generateSequenceGCContent(vals[0],.6))
    sequences.append(ancestor)

    #once the generation parameter runs out of steam, break from the while loop
    while vals[1] > 0:
        if possibleDuplication():
            vals[1] = vals[1] - 1
            sequences.append(possibleMutationDuringDuplication(ancestor.sequence))
        else:
            for sequence in sequences:
                if possibleMutation():
                    sequence.randomModify()

    return


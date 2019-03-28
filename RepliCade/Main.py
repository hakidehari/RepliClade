"""

This is the Main file for Repliclade where the program will be run
All dependencies and classes will be compiled here and ran accordingly

"""

import random
from Generator import Generator
from Sequence import Sequence
from Connector import Connector
from ParseFasta import ParseFasta
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import os


'''
    Performs multiple sequence alignment after simulation is run
    using clustalW2, one of the fastest MSA algorithms out there.
    Once finished, displays the sequences and their alignment
'''
def performMSA():
    clustalw_exe = r"/applications/clustalw-2.1-macosx/clustalw2"
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile="influenza.fasta")
    stdout, stderr = clustalw_cline()
    align = AlignIO.read("influenza.aln", "clustal")
    print(align)


'''Gathers sequences from GenBank'''
def getSequences():
    con = Connector()
    seq = con.getGeneData()
    return seq


'''returns the user specified DNA length and simulation time (in iterations)'''
def simulationParameters(DNALength,simulationTime):
    return [DNALength, simulationTime]


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
    if r < 0.2:
        sequence = Sequence(sequence.randomModify())
    return sequence


'''will run through a probability model and return whether the sequence should mutate or not'''
def possibleMutation():
    r = random.random()
    if r < .05:
        return True
    return False

'''Prints the sequences returned from the simulation'''
def printSequences(sequences):
    ar = [];
    for sequence in sequences:
        print(sequence.sequence)
        ar.append(len(sequence.sequence))
    print(len(sequences))
    print(ar)

'''
    Run simulation given two user specified parameters returned from the simulationParameters() function
    Currently runs smoothly with 10k iterations but slows down tremendously after 50k iterations
    What causes this is the fact that too many sequences are being created and the simulation cannot
    keep up with the iteration through the sequences.  I have brainstormed possible solutions.
    
    We hope to implement MSA soon and construct a phylogenetic tree of the sequences simulated.
'''
def runSimulationRandom():
    #vals[0] is the DNA length and vals[1] is the simulation time in (units here)
    vals = simulationParameters(500, 10000)
    seqLen = vals[0]
    runTime = vals[1]
    sequences = []
    g = Generator()
    a = g.generateSequenceGCContent(seqLen, 0.6)
    ancestor = Sequence(a)
    sequences.append(ancestor)
    #use a previous array of sequences to prevent throttling as much as possible
    prevIteration = sequences

    #once the generation parameter runs out of steam, break from the while loop
    while runTime > 0:
        for i in range(0, len(prevIteration)):
            #check if duplication should happen
            if possibleDuplication():
                #append with a possible error in duplication
                sequences.append(possibleMutationDuringDuplication(sequences[i]))
            #check if indel insertion should happen
            if possibleIndelInsertion():
                #insert an indel somewhere into the sequence
                sequences[i] = Sequence(sequences[i].insertIndel())
            else:
                #else check the probability of a nucleotide mutation happening
                if possibleMutation():
                    #if so, modify it and save it to the sequence list
                    sequences[i] = Sequence(sequences[i].randomModify())

        prevIteration = sequences
        runTime -= 1

    printSequences(sequences)
    return sequences


'''
    Run Simulation off of 8 influenza genomes (this can be anything in the future, we used influenza
    for testing purposes)
    runTime is set to 10 thousand for now.  The proper probability models still need to be implemented
    Mutation and duplication is still happening more often than we would hope.
'''
def runSimulationGenome():
    influenza = getSequences()
    runTime = 10000
    sequences = influenza
    prevIteration = influenza

    while runTime > 0:
        for i in range(0, len(prevIteration)):
            #check if duplication should happen
            if possibleDuplication():
                #append with a possible error in duplication
                sequences.append(possibleMutationDuringDuplication(sequences[i]))
            #check if indel insertion should happen
            if possibleIndelInsertion():
                #insert an indel somewhere into the sequence
                sequences[i] = Sequence(sequences[i].insertIndel())
            else:
                #else check the probability of a nucleotide mutation happening
                if possibleMutation():
                    #if so, modify it and save it to the sequence list
                    sequences[i] = Sequence(sequences[i].randomModify())

        prevIteration = sequences
        runTime -= 1

    printSequences(sequences)
    parseObj = ParseFasta()
    parseObj.writeToFasta(sequences)

    performMSA()

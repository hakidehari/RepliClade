"""

This is the Main file for Repliclade where the program will be run
All dependencies and classes will be compiled here and ran accordingly

"""

import random
from Generator import Generator
from Sequence import Sequence
from Connector import Connector
from EyeOh import EyeOh
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import Phylo
from ProbabilityModels import *


'''
    Performs multiple sequence alignment after simulation is run
    using clustalW2, one of the fastest MSA algorithms out there.
    Once finished, displays the sequences and their alignment
'''
def performMSA():
    print("Performing Multiple Sequence Alignment.......")
    clustalw_exe = r"/applications/clustalw-2.1-macosx/clustalw2"
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile="influenza.fasta")
    stdout, stderr = clustalw_cline()
    align = AlignIO.read("influenza.aln", "clustal")
    print(align)
    print("Alignment Complete.")


'''
    Creates and displays a phylogenetic tree based on the Multiple Sequence Alignment scores
'''
def displayTree():
    tree = Phylo.read("influenza.dnd", "newick")
    tree.rooted = False
    Phylo.draw(tree)


'''Gathers sequences from GenBank'''
def getSequences(gene):
    con = Connector()
    seq = con.getGeneData(gene)
    return seq


'''returns the user specified DNA length and simulation time (in iterations)'''
def simulationParameters(DNALength,simulationTime):
    return [DNALength, simulationTime]


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
    influenza = getSequences(None)
    print("Simulating..........")
    runTime = 50000
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

    parseObj = EyeOh()
    parseObj.writeToFasta(sequences)
    print("Simulation Complete.")

    performMSA()
    displayTree()


#find percentage of similarity between 2 sequences
def alikeness(s1, s2):
    count = 0
    for i in range(0, len(s1)):
        if s1[i] == s2[i]:
            count += 1
    return count / len(s1)


def sequencify(arr):
    for i in range(0, len(arr)):
        arr[i] = Sequence(arr[i])
    return arr


#clean duplicate sequences in sequence list
def cleanSequences(seqArray):
    returnArr = []
    returnArr.append(seqArray[0].sequence)
    for i in range(1, len(seqArray)):
        if seqArray[i].sequence not in returnArr:
            returnArr.append(seqArray[i].sequence)
    returnArr = sequencify(returnArr)
    return returnArr



'''
    Simulation starting with a single ancestor.  Will be trying to implement some sort of extinction
    mechanism
'''
def runSimulationSingleAncestor():
    ancestor = getSequences("KT388711")[0]
    print("Simulating")
    #run time in generations
    runtime = 15
    #the rate of reproduction for the influenza A virus
    r0 = 1.5
    sequences = [ancestor]
    current = sequences
    newCurrent = []
    while runtime > 0:
        for i in range(0, len(current)):
            newSequence = influenzaMutate(current[i])
            sequences.append(newSequence)
            newCurrent.append(newSequence)
            if r0 == 3:
                newSequence = influenzaMutate(current[i])
                sequences.append(newSequence)
                newCurrent.append(newSequence)
        if r0 == 1.5:
            r0 = 3
        else:
            r0 = 1.5
        current = newCurrent
        newCurrent = []
        runtime -= 1

    sequences = cleanSequences(sequences)
    parseObj = EyeOh()
    parseObj.writeToFasta(sequences)
    print("Simulation Complete.")

    performMSA()
    displayTree()














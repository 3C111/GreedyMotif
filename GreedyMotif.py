#!/usr/bin/python

## The algorithm iteratively finds k-mers in each string of DNA
## After i-1 kmers are found in the first i-1 strings of DNA a profile matrix is constructed
## and the Profile most probable kmer is selected from the i-th string

import numpy as np
from Profilemostprobablekmer import mostprobablekmer

def NormalizeArray(array, n):
    # Normalize a counts array to n

    for i in range(0, len(array)):
        if array[i] != 0:
            array[i] = array[i]/float(n)
        else:
           array[i] = 0 

    return array

def buildProfileMatrix(Motifs):
    noRows = len(Motifs)
    noColumns = len(Motifs[0])

    # Initalise the counts arrays
    countA = [0]*noColumns
    countC = [0]*noColumns
    countG = [0]*noColumns
    countT = [0]*noColumns
    

    # Create the counts matrix
    for motif in Motifs:
        for i in range(0, noColumns):
            if motif[i] == 'A':
                countA[i] +=1
            if motif[i] == 'C':
                countC[i] +=1
            if motif[i] == 'G':
                countG[i] +=1
            if motif[i] == 'T':
                countT[i] +=1

    # Normalize the counts to get the profile matrix using the helper function NormalizeArray

    profilematrix = []
    profilematrix.append(NormalizeArray(countA, noRows))
    profilematrix.append(NormalizeArray(countC, noRows))
    profilematrix.append(NormalizeArray(countG, noRows))
    profilematrix.append(NormalizeArray(countT, noRows))

    return profilematrix

def Score(Motifs):
   # The score finds the most common nucleotide in each column and counts the amount of nucleotides 
   # that are not the most common one

    noRows = len(Motifs)       # Number of motifs
    noColumns = len(Motifs[0]) # Length of the individual Motif
    
    # Determine the dominant nucleotide in each column
    
    dominant_nucleotide = []

    for i in range(0, noColumns):    
        nucleo_occurrence = [0, 0, 0, 0]  # A, C, G, T
        
        for motif in Motifs: 
            if motif[i] == 'A':
                nucleo_occurrence[0]+=1
            if motif[i] == 'C':
                nucleo_occurrence[1]+=1
            if motif[i] == 'G':
                nucleo_occurrence[2] +=1
            if motif[i] == 'T':
                nucleo_occurrence[3] +=1
        # Add to dominant nucleotide the most common peptide    
        dominant_nucleotide.append(np.argmax(nucleo_occurrence))
    
    # Count the amount of nucleotides that are not the dominant nucleotide
    nucleo = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

    score = 0
    for j in range(0, noColumns):
        dominant_nucleo = nucleo[dominant_nucleotide[j]]
        for motif in Motifs:
            if motif[j] is not dominant_nucleo:
                score+=1
       
    return score

def Greedymotifsearch(Dna, k, t):

    # t is the amount of strands in Dna
    # Create the inital bestmotif matrix from the first kmers in each strand Dna

    BestMotifs = []

    for i in range(0, len(Dna)):
        strand = Dna[i]
        BestMotifs.append(strand[0:k])

    firstStrand = Dna[0]

    for j in range(0, len(firstStrand)-k+1):
        # Start with a kmer in the first strand
        Motifs=[]
        Motifs.append(firstStrand[j:j+k])
         
        for i in range(1,t):
            strand = Dna[i]
            # Create a profile with the motifs in Motifs 
            profile = buildProfileMatrix(Motifs)
            # Add the most probable motif based on the profile
            Motifs.append(mostprobablekmer(strand, k, profile))

        if  Score(Motifs) <= Score(BestMotifs):
            BestMotifs = Motifs
            print Score(BestMotifs)
    return BestMotifs

# Sample Input: k = 3, t = 5
# DNA= ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']

# Sample Output: CAG, CAG, CAA, CAA, CAA
k = 3
t = 5

Dna=['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']


motifs= Greedymotifsearch(Dna, k, t)
 
for motif in motifs:
    print  motif

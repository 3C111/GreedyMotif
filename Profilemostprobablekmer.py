#!/usr/bin/python 

## Find the most probable k-mer in Text from the profile matrix
## Calculate the probability of every k-mer in Text and select the one with
## the highest probability

def Probability(Pattern, Profile):

    # Divide the Profile in lists per nucleotide
    A = Profile[0]
    C = Profile[1]
    G = Profile[2]
    T = Profile[3]
 
    # Initialise parameters
    probability = 1.0
   
    # Loop over Pattern multiply corresponding probabilities
    for i in range(0, len(Pattern)):
        if Pattern[i] == 'A':
            probability*=A[i]
        if Pattern[i] == 'C':
            probability*= C[i]
        if Pattern[i] == 'G':
            probability*=G[i]
        if Pattern[i] == 'T':
            probability*=T[i]

    
    return probability 


def mostprobablekmer(Text, k, Profile):
    # Initialise parameters
    kmer = ''
    probability = 0

    # Look over Text to test each pattern
    for i in range(0, len(Text)-k+1):
        pattern = Text[i:i+k]
         
        Probability_pattern = Probability(pattern, Profile)
         
        if probability < Probability_pattern:
            probability = Probability_pattern             
            kmer = pattern
    if kmer == '':
        kmer = Text[0:k]        
    return kmer

# Sample Input Text = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
# k = 5
# Profile = [[ 0.2, 0.2, 0.3, 0.2, 0.3], [0.4, 0.3, 0.1, 0.5, 0.1], [0.3, 0.3, 0.5, 0.2, 0.4] [0.1, 0.2, 0.1, 0.1, 0.2]]

# Sample Output: CCGAG



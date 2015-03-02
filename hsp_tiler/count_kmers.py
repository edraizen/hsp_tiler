#!/usr/bin/env python2.7
#Author: Eli Draizen
#BME 205 Assignment 4
#Date: 10/31/12

"""Module to count kmers in a given FASTA files

Adapted from Prof. Karplus' Markov Assignment for BME 205 f12
"""

from collections import Counter
from read_fasta import read_fasta

IUPAC_AA = "ACDEFGHIKLMNOPQRSTVWY"
IUPAC_N = "ACTG"

def count_kmers(*args, **kwds):
    #inFile, alphabet, [m,] k, step=1, stopCharacter="="
    """Count the frequency of k-mers in a given file and alphabet.
    -inFile, file-like object containing FASTA sequences
    -alphabet, the alphabet the FASTA sequence uses. Cahracters no in the 
        alphabet are ignored,
    -m, lower bound for kmer length, optional
    -k, upper bound for kmer length, inclusive
    -step, step size. Default 1
    -stopCharacter, define the stop character used signify the start and end of 
        sequences. Optional.
    -saveLength, saves legnth of each sequence and return the sum

    Kwds-
    normalize=False
    """
    if len(args) < 3:
        raise RuntimeError("count_kmers requires at least 3 arguments")
    if len(args) == 4:
        m = args[2]
        k = args[3]
    else:
        m = k = args[2]
        
    if m and m<1:
        raise RuntimeError("Cannot generate a kmer less than 1")
    if m>k:
        raise RuntimeError("The lower bound must be lower than or equal to the \
upper bound.")

    inFile = args[0]
    alphabet = args[1]
    step = kwds.get("step", 1)
    stopCharacter = kwds.get("stopCharacter", "=")
    saveLength = kwds.get("saveLength", False)
    
    if type(inFile) == type([]) and inFile[0].__class__.__name__ == "Sequence":
        #inFile is already parsed FASTA
        seqs = inFile
    else:
        #parse incoming fasta
        seqs = read_fasta(inFile, alphabet)
    
    counts = Counter()
        # counts['ABC'] is the number of times the 3-mer 'abc' occurs 
        # in the input.
        # Counter used to make counts of unseen k-mers be 0
        
    if saveLength:
        #length of each sequence
        numSeqs = []
    else:
        #number of all FASTA sequences
        numSeqs = 0
    
    for fastaSeq in seqs:
        if saveLength:
            numSeqs.append(len(fastaSeq.sequence))
        else:
            numSeqs += 1
        
        for palindromeLength in range(m, k+1, step):
            #Pad sequence with non seq character
            if palindromeLength>1:
                start = stopCharacter*(palindromeLength-1)
                end = stopCharacter*(palindromeLength-1)
            else:
                start = ""
                end = stopCharacter
            seq = "{}{}{}".format(start, fastaSeq.sequence, end)
            for i in range(len(seq)-palindromeLength+1):
                kmer = seq[i:i+palindromeLength]
                counts[kmer] += 1

    if kwds.get("normalize", False):
        total = sum(counts.values())
        for kmer, count in counts.items():
            counts[kmer] = counts[kmer]/total
        
    return counts, numSeqs

def read_kmers(self, kmerFile):
    """Read in a file that contains kmers and their frequencies. The length of 
kmer is figure out. File must have the kmer and count spearated by a single space.
    
    -kmerFile, file-like object with kmer counts
    """
    
    #Initialize a kmer couts counter, sets unseen kmers to 0
    kmer_counts = Counter()
    
    #Determine length of kmer
    k = None
    
    for line in kmerFile:
        line = line.rstrip()
        
        fields = line.split()
        
        #Line must have the kmer and count spearated by a single space.
        if len(fields) > 2:
            raise RuntimeError("Invalid kmerFile! There is more the one field on \
line {}".format(line))
            
        #The kmer is the first field
        kmer = fields[0]
       
        #The second field must be the count and must be the same length as other kmers
        if k and not len(kmer) == k:
            raise RuntimeError("Invalid kmerFile! The second field must be the count \
and must be the same length as other kmerson line {}".format(line))
        else:
             k = len(kmer)
             
        #The second field must be an integer
        if fields[1].isdigit():
            count = int(fields[1])
        else:
            raise RuntimeError("Invalid kmerFile! Count must bean int on \
line {}".format(line))
            
        kmer_counts[kmer] = count
    
    return kmer_counts

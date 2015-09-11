#!/usr/bin/env python2.7
#Author: Eli Draizen
#BME 205 Assignment 3
#Date: 10/23/12

"""Object to intialize and train a k-order Markov Model
by counting kmers in a given sequence

Use by calling:
import Markov
markov_model = Markov.Markov(k, alphabet, countFile)
"""

import re
import math
from collections import defaultdict, Counter
from itertools import product, groupby

class Markov(object):
    def __init__(self, kmer_counts, k, alphabet, stopCharacter="=", pseudocount=1.0, ):
        """Initialize a Markov object to train and score a k-order markov models. 
        The markov model is trained by passing in a dictionary with kmers as keys and 
        their frequency as values, an alphabet, and length of kmer.

        Arguments:
        k:Markov order. Split a sequence into strings of length k.
        alphabet: string, characters used in sequence.
        kmer_counts: dictionary with kmers as keys and their frequency as values,
        stopCharacter: Append a certain chachter to end of strings to signigy stop. Defualt = '='.
        kmerFile: file-like object containing a pre-built log probability table. Optional.
        psuedocount: the psudocount to normalize data that hasn't been seen before,
        """
        #Save parameters
        self.k = k
        self.alphabet = alphabet
        self.stopCharacter = stopCharacter

        #Train Markov model

        #Define log probability table:
        #Create the dictionary of keys as kmers and values as the probability
        #Unseen kmers are set to 0 at first, but corrected with pseudocounts
        self.log_prob_table = defaultdict(float)
        
        #Create regex that will find stop charters in the middle of strings
        pat = re.compile("[{}]+{}+[{}]+".format(re.escape(alphabet),
                                                re.escape(stopCharacter),
                                                re.escape(alphabet)))
        #Create every possible permutation of the alphabet with length k
        for kmer in product(alphabet+stopCharacter, repeat=k):
            #Turn kmer back into a string instead of a tuple
            kmer = "".join(kmer)
            
            #See if there is a stop in middle followed by something other than stop, 
            #if so skip
            if pat.search(kmer):
                continue
                
            #Add the kmer to the log table with the pseudocount and the count, if seen
            self.log_prob_table[kmer] = pseudocount + kmer_counts[kmer]
            
        #Normalize the counts to sum to 1 over each "condition".
        #Group the kmers together by condition
        for condition, kmers in groupby(sorted(self.log_prob_table.keys()), key=lambda x:x[:-1]):
            #Sum up all kmers with same condition
            total = 0
            
            kmers = list(kmers)
            for kmer in kmers:
                total += self.log_prob_table[kmer]
            
            #Normalize by dividing count by total
            for kmer in kmers:
                prob = float(self.log_prob_table[kmer])/float(total)
                self.log_prob_table[kmer] = -math.log(prob, 2)

    def score(self, sequence):
        """Score a given sequence based on the trained markov model
        Input:
        sequence: string that contains characters witht he same alphabet
        cost: boolean, divide by sequence length?
        """
        #Score of the sequence (sum of encoding of each character)
        sequence_score = 0
        
        #Pad sequence with non seq character
        if self.k>1:
            start = self.stopCharacter*self.k
            end = self.stopCharacter*self.k
        else:
            start = ""
            end = self.stopCharacter
        seq = "{}{}{}".format(start, sequence, end)
        
        for i in range(len(sequence)-self.k+1):
            kmer = sequence[i:i+self.k]
            character_score = self.log_prob_table[kmer]
            #print character_encoding_cost
            sequence_score += character_score

        return sequence_score/float(len(sequence))




    
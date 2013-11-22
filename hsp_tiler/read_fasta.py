#Author: Eli Draizen
#File: read_fasta
#Date: 10/10/12

import sys, re
import string

class Sequence(object):
    def __init__(self):
        #Initialize an empty sequence
        self.name = ""
        self.description = ""
        self.sequence = ""

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        space = len(self.description)>0 and " " or ""
        tmp = ">{}{}{}\n".format(self.name, space, self.description)
        for i in range(0, len(self.sequence), 60):
            tmp += "{}\n".format(self.sequence[i:i+60])
        return tmp[:-1]

def read_fasta(inFile, alphabet=None):
    """Read in fasta file and return sequence objects with the fasta data
    
    inFile - a file-like object contianing the fasta data
    """
    #Set up sequence and flags for sequence and quality
    fastaSeq = None
    empty = True
    
    #Create alphabet for regex based on IUPAC amino acid and nucleic acids codes
    if alphabet:
        pattern = "[{}]+".format(re.escape(alphabet))
        alphabet = re.compile(pattern, re.IGNORECASE)
    else:
        alphabet = re.compile(r"[A-Z]+", re.IGNORECASE)
        delete = ""
    
    for i, line in enumerate(inFile):
        line = line.rstrip()
        if fastaSeq and line.startswith(">"):
            yield fastaSeq
        
        if line.startswith(">"):
            #Start new sequence and get id and comment
            fastaSeq = Sequence()
            
            #Get name and desciption
            try:
                fastaSeq.name, fastaSeq.description = line[1:].split(None, 1)
            except:
                fastaSeq.name = line[1:]

                
            empty = False
        elif fastaSeq is None:
            #Invalid FASTA, terminate program
            print >> sys.stderr, "Error: Fasta file must begin with '>' \
on line '{}'".format(line)
            sys.exit(1)
        else:
            #Append line to sequence
            #Translate may be faster, but I kept on removing the wrong characters
            for seq in alphabet.findall(line):
                fastaSeq.sequence += seq.upper()
    
    #Return last fasta entry if exists            
    if fastaSeq:
        yield fastaSeq
    
    #Warn if empty or has errors
    if empty:
        print >> sys.stderr, "WARNING: fastaFile has no sequences."

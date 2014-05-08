#Author: Eli Draizen
#File: read_fasta
#Date: 10/10/12

#Standard libraries
import sys, re
import string

#Custom libraries
from Sequence import DNASequence 

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
            fastaSeq = DNASequence()
            
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
            for seq in alphabet.findall(line):
                fastaSeq.sequence += seq.upper()
    
    #Return last fasta entry if exists            
    if fastaSeq:
        yield fastaSeq
    
    #Warn if empty or has errors
    if empty:
        print >> sys.stderr, "WARNING: fastaFile has no sequences."

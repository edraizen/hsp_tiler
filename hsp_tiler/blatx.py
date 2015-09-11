#!/usr/local/bin/python
# Author: Eli Draizen
# Date: 16-3-2014
# File: blatx.py

#Standard Libraries
import argparse
import sys

#Required Libraries
import subprocess

#Custom Libraries
from read_fasta import read_fasta
from Sequence import Sequence

def blatx(dna, protein, outfile, run=False):
	sixframeFileName = "sixframe.blat.{}".format(dna.name)
	with open(sixframeFileName, "w") as sixframeFile:
		for fasta in read_fasta(dna):
			for frame, protein in fasta.sixframe():
				protein = protein.replace("_", "*")
				protSeq = Sequence(name="{}_f{}".format(fasta.name, frame),
					               sequence=protein)
				print >> sixframeFile, protSeq

	if run:
		p = Popen(["blat", "-t=prot", "-q=prot", "-out=pslx", sixframeFileName, protein], 
			      stdout=PIPE, stdin=PIPE, stderr=PIPE)
		p.communicate()


def parse_args():
    """Parsing command line options
    """
    parser = argparse.ArgumentParser(description="")
    
    #Define Input
    parser.add_argument("-f", "--fasta", 
                        required=True, 
                        type=argparse.FileType('r'),
                        help="Fasta file containing sequences")
    parser.add_argument("-p", "--protein", 
                        help="File containing protein sequences")
    parser.add_argument("--run", 
                        default=False,
                        action="store_true",
                        help="Run BLAT!")
 
    #Define output
    parser.add_argument("-o", "--outfile",
    					default="/dev/null",
                        required=False,
                        help="File to save corrected sequences")

    #Parse args
    return parser.parse_args()

if __name__ == "__main__":
    #Parse args
    args = parse_args()

    blatx(args.fasta, args.protein, args.outfile)




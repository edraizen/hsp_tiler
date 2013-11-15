#Standard Libraries
import sys
import re
import argparse
import csv
from itertools import izip

#Required Libraries
from Bio.Blast.Applications import NcbiblastxCommandline
import prettyplotlib as ppl
import numpy as np
from prettyplotlib import plt

#Custom imports
from read_fasta import read_fasta
import count_kmers
from read_codon_usage import read_codon_frequencies

"""Score the output from hsp_tiler and optionally create pretty graphics"""

def score_blastx(originalFilename, updatedFilename):
	"""Compare the blastx bitscores of the original fasta file and hsp-tiler
	Parameters;
	originalFilename - filename of original fasta files used to run hsp_tiler
	updatedFilename - filename of hsp_tiler output
	"""
	#Calculate original scores
	original_scores = get_blastx_scores(originalFilename)

	#Calculate updated scores
	updated_scores = get_blastx_scores(updatedFilename)

	return original_scores, updated_scores

def write_scores(fasta, original_scores, updated_scores, outfile):
	print >> outfile, "Contig\tOld Score\tNew Score\thsp_tiler_score"
	for sequence, oldScore, newScore in izip(read_fasta(fasta), original_scores, updated_scores):
		score = newScore/oldScore
		print >> outfile, "{}\t{}\t{}\t{}".format(sequence.name, oldScore, newScore, score)

def read_scores(score_file):
	"""Read in a score file written by the write_scores functions. Useful for re running script
	to edit graphs or view the in iPython notbook

	Parameters:
	score_file, a file-like object containing tab delimeted score informations

	Output:
	Three lists, old scores, new scores, and the hsp-tile, score
	"""
	old_scores = []
	new_scores = []
	scores = []
	for row in csv.DictReader(score_file, dialect="excel-tab"):
		old_scores.append(float(row["Old Score"]))
		new_scores.append(float(row["New Score"]))
		scores.append(float(row["hsp_tiler_score"]))
	return old_scores, new_scores, scores

def markov_score(original, updated, method=0):
	if method == 0:
		#Create log probability table
		counts, num_seqs = count_kmers.count_kmers(original, count_kmers.IUPAC_N, 3)
	else:
		#Codon usage for C. elegans
		counts = read_codon_frequencies(6239)

	#Train 3rd order markov model, based on codons 
	markov_model = Markov.Markov(kmer_count, 3, count_kmers.IUPAC_N)

    #original_score = markov_model.score(contig.sequence)
    #original_scores.append(original_score)
    #print "{} score: {}".format(queryID, original_score)

def lev_distance(s1, s2):
	"""Copied from ... just to test, if works will correctly cite"""
	if len(s1) > len(s2):
		s1,s2 = s2,s1
	distances = range(len(s1) + 1)
	for index2,char2 in enumerate(s2):
		newDistances = [index2+1]
		for index1,char1 in enumerate(s1):
			if char1 == char2:
				newDistances.append(distances[index1])
			else:
				newDistances.append(1 + min((distances[index1],
									distances[index1+1],
									newDistances[-1])))
		distances = newDistances
	return distances[-1]

	

def get_blastx_score(seq):
	#Save contig to be used by blastx
	blastx_seq = open("blastx_seq.fa", "w+")
	blastx_seq.write(str(seq))
	blastx_seq.close()

	#Run blastx
	cmd = "{} -db {} -query blastx_seq.fa -max_target_seqs 1 -outfmt '6 bitscore'".format(blastx, blastdb)
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
	out, err = process.communicate()

	#Cleanup
	os.remove("blastx_seq.fa")

	return int(out.strip())

def get_blastx_scores(fastaFile):
	blastFile = "{}.blast.txt".format(fastaFile)
	print >> sys.stderr, blastFile
	try:
		blastFile = open(blastFile)
	except:
		blastx_cline = NcbiblastxCommandline(query=fastaFile, db=blastdb, max_target_seqs=1,
	                                         outfmt=6, num_threads=11, out=blastFile)
		stdout, stderr = blastx_cline()
		blastFile = open(blastFile)

	import csv
	scores = [float(field[11]) for field in csv.reader(blastFile, delimiter="\t")]
	#{field[0]: float(field[11]) for field in csv.reader(blastFile, delimiter="\t")}

	return scores
   
def scatter(original, updated, xlab=None, ylab=None, main=None):
	fig, ax = plt.subplots()
	ax.set_title(main)
	ax.set_xlabel(xlab)
	ax.set_ylabel(ylab)
	ppl.scatter(ax, original, updated)
	plt.show()

def histogram(original, updated):
	fig, ax = plt.subplots()

	ohist, obins = np.histogram(original, bins=50)
	width = 0.7*(obins[1]-obins[0])
	center = (obins[:-1]+obins[1:])/2
	ppl.bar(ax, center, ohist, align='center', width=width)

	uhist, ubins = np.histogram(updated, bins=50)
	width = 0.7*(ubins[1]-ubins[0])
	center = (ubins[:-1]+ubins[1:])/2
	ppl.bar(ax, center, uhist, color='r', align='center', width=width)

	plt.show()

def parse_args():
	parser = argparse.ArgumentParser(description="Takes a fasta file of sequences and a BLASTX annotation of that file in xml format.  Attempts to tile Hsps for the highest scoring hit for each sequence, correcting frameshifts in order to improve subsequent annotations.")
	# name of fasta file 
	parser.add_argument("-o", "--original", 
	                    required=True, 
	                    type=argparse.FileType('r'),
	                    help="Fasta file containing sequences")
	# name of blast output in tab delimated (standard)
	parser.add_argument("-u", "--updated", 
	                    required=True, 
	                    type=argparse.FileType('r'),
	                    help="Blastx xml file containing annotations for sequences")
	#blast for evalue comparison
	parser.add_argument("-b", "--blastdb",
	                    required=True,
	                    help="Location of blast db to caluate evalue changes")
	parser.add_argument("--blastx",
	                    required=False,
	                    default="blastx",
	                    help="Location of blastx executible, if not in users PATH")
	#Define output
	parser.add_argument("-s", "--outfile",
	                    type=argparse.FileType('wt'),
	                    default=sys.stdout,
	                    help="File to save corrected sequences")
	return parser.parse_args()

if __name__ == "__main__":
	args = parse_args()
	blastdb = args.blastdb
	blastx = args.blastx
	original_scores, updated_scores = score_blastx(args.original.name, args.updated.name)
	write_scores(args.updated, original_scores, updated_scores, args.outfile)
	scatter(original_scores, updated_scores)
#Major Bug: Script does not account for filters set up by hsp-tiler
#Fix by using the score that started the tile... 


#Standard Libraries
import sys
import re
import argparse
import csv
import subprocess
from itertools import izip, tee
from collections import defaultdict

#Required Libraries
from Bio.Blast.Applications import NcbiblastxCommandline
#import prettyplotlib as ppl
#import numpy as np
#from prettyplotlib import plt

#Custom imports
from read_fasta import read_fasta
import count_kmers
from read_codon_usage import read_codon_frequencies

"""Score the output from hsp_tiler and optionally create pretty graphics"""

def get_original_scores(fastaFile, saveGI=False):
	#Parse out gi number and score from each updated contig
	seqPattern = re.compile("\[GI\=(.*?)\;S\=(.*?)\]")
	originalScores = []
	not_corrected = []
	giFileName = "{}.GILIST.txt".format(fastaFile)
	ids = {}
	with open(fastaFile) as f, open(giFileName, "w") as giFile:
		for seq in read_fasta(f):
			seqInfo = seqPattern.match(seq.description)
			if not seqInfo:
				raise RuntimeError("Sequnces must contain HSP-Tiler headers.")

			gi = seqInfo.group(1)
			if saveGI and not gi == "None":
				print >> giFile, gi 
			ids[f.name] = gi
			originalScores.append(float(seqInfo.group(2)))

	if saveGI:
		return originalScores, giFileName, ids
	return originalScores, ids

def score_blastx(fastaFile, blastdb, makegilist=False, blastpath=""):
	"""Calculate updated HSP-Tiler scores. The sequences that were used to 
	create the tile are used to filter the blast database in order to remove 
	hits that were filtered during HSP-Tiler.

	Input:
	fastaFile - path to HSP-Tiler output FASTA
	blastdb - path to the BLAST db
	makegilist - bool. Create a new databse filtered by sequenced used to start the tile
	blastx - path to BLASTX executable
	"""
	
	if makegilist:
		if blastpath and not blastpath.endswith("/"):
			blastpath += "/"
		newblastdb = "{}.gi.blastdb".format(fastaFile)
		cmd = "{}blastdb_aliastool -dbtype prot -gilist {} -db {} -out {}".format(blastpath, giFileName, blastdb, newblastdb)
		process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
		out, err = process.communicate()
		blastdb = newblastdb

	print "Running BLAST"
	blastFile = "{}.blast.txt".format(fastaFile)
	blastx_cline = NcbiblastxCommandline(query=fastaFile, db=blastdb, max_target_seqs=10,
	                                     outfmt=6, num_threads=4, out=blastFile)
	print >> sys.stderr, "Running BLASTX"
	stdout, stderr = blastx_cline()

	with open(blastFile) as b:
		updated_scores = read_blastx_bitscores(b)

	return updated_scores

def read_blastx_bitscores(resultsFile, ids):
	#Process blastx ouput
	#return [float(field[11]) for field in csv.reader(resultsFile, delimiter="\t")]
	scores = []
	previousID = None
	for field in csv.reader(resultsFile, delimiter="\t"):
		if field[0] == previousID and field[1] == ids[field[0]]:
			scores.append(float(field[11]))
		previousID = field[0]
	return scores

def write_scores(fasta, original_scores, updated_scores, outfile):
	"""Write scores to File

	Input:
	fasta - file-like object containing FASTA sequences
	original_scores - a list of scores
	updated_scores - a list of scores
	"""
	print >> outfile, "Contig\tOld Score\tNew Score\thsp_tiler_score"
	for sequence, oldScore, newScore in izip(read_fasta(fasta), original_scores, updated_scores):
		if oldScore != 0.0:
			score = newScore/oldScore
		else:
			score = 0.0
		print >> outfile, "{}\t{}\t{}\t{}".format(sequence.name, oldScore, newScore, score)

def read_scores(score_file):
	"""Read in a score file written by the write_scores function. Useful for re running script
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
   
def scatter(original, updated, xlab=None, ylab=None, main=None):
	fig, ax = plt.subplots()
	ax.set_title(main)
	ax.set_xlabel(xlab)
	ax.set_ylabel(ylab)
	ppl.scatter(ax, original, updated)
	plt.show()

def histogram(original, updated, bins=None):
	"""Plot a histogram of score improvements (updated-origianl)

	Input:
	original - list of original scores
	updated - list of updates scores in same order as original
	bins - number of bins to represent improvements
	"""
	import numpy as np
	#Lengths of score lists must be identical, assume in same order
	assert len(original) == len(original)

	#Set up bins:
	if bins is None or bins != 0:
		imoprovements = {(-1,-1):0}
		for i in range(0, len(original), bins):
			improvements[(0,i+bins)] = 0
	else:
		improvements = {(-1,-1):0, (-5,0):0, (0,1):0, (1,25):0, (25,50):0, (50,75):0, (75,100):0, (100,125):0, (125,150):0, (150,200):0, (200,300):0, (300,400):0, (500,10000):0} #defaultdict(int)
	
	#Calcualte improvements
	for o, u in izip(original, updated):
		if o>u: 
			improvements[(-1,-1)] += 1
			continue
		for lower, upper in improvements:
			if lower <= int(u-o) < upper:
				improvements[(lower,upper)] += 1
				break
	keys = sorted(improvements.keys(), key=lambda x:x[0])
	values = [improvements[r] for r in keys]

	fig, ax = plt.subplots()
	#ppl.hist(ax, improvements.values())

	width = 1.0
	ax.set_xticks(np.arange(len(improvements)))
	ax.set_xticklabels([l for l, u in keys])
	ppl.bar(ax, np.arange(len(improvements)), values, align='center', width=width, log=True)

	plt.show()

def parse_args():
	parser = argparse.ArgumentParser(description="Score HSP-Tiler output FASTA files")
	parser.add_argument("-f", "--fasta", 
	                    required=True,
	                    type=argparse.FileType('r'),
	                    help="Blastx xml file containing annotations for sequences")
	
	subparsers = parser.add_subparsers(help='Run BLAST, load BLAST RESULTS, or load score file')

	# Run BLASTX
	runblast = subparsers.add_parser('run', help='Run/preprocess BLAST')
	runblast.add_argument("-b", "--blastdb",
	                    required=True,
	                    help="Location of blast db to calcuate evalue changes")
	runblast.add_argument("--blastpath",
	                    required=False,
	                    default="",
	                    help="Location of blast executibles, if not in user's PATH")
	runblast.add_argument("-g", "--makegilist",
	                    required=False,
	                    action="store_true",
	                    help="Score sequenes from filtered blastdb based on gi number")


	# Read BLASTX reults
	readblast = subparsers.add_parser('results', help='Read BLASTX results')
	readblast.add_argument("-r", "--results",
					       type=argparse.FileType('r'),
					       help="BLAST output file with updated scores")

	#Define output
	parser.add_argument("-o", "--outfile",
	                    type=argparse.FileType('wt'),
	                    default=sys.stdout,
	                    help="File to save corrected sequences")
	return parser.parse_args()

if __name__ == "__main__":
	#Parse arguments
	args = parse_args()

	if hasattr(args, "results"):
		original_scores, ids = get_original_scores(args.fasta.name)
		updated_scores = read_blastx_bitscores(args.results, ids)
	elif hasattr(args, "run") and args.makegilist:
		original_scores, giFileName = get_original_scores(fastaFile, saveGI=True)
		updated_scores = score_blastx(args.fasta.name, args.blastdb, makegilist=giFileName, blastpath=args.blastpath)
	else:
		original_scores, ids = get_original_scores(args.fasta.name)
		updated_scores = read_blastx_bitscores(args.fasta.name, ids)

	#Save scores
	write_scores(args.fasta, original_scores, updated_scores, args.outfile)

	#Plot
	#scatter(original_scores, updated_scores)

#!/usr/bin/python2.7
#Author: Eli Draizen 
#File: hsp_score.py

#Standard Libraries
import sys
import re
import argparse
import csv
import subprocess
from itertools import izip, tee
from collections import defaultdict

#Required Libraries
try:
	from Bio.Blast.Applications import NcbiblastxCommandline
	hasBLAST = True
except ImportError:
	hasBlast = False

#Custom imports
from read_fasta import read_fasta

"""Score the output from hsp_tiler"""

def get_original_scores(fastaFile):
	#Parse out gi number and score from each updated contig
	seqPattern = re.compile("\[GI\=(.*?)\;S\=(.*?)\]")
	originalScores = []
	not_corrected = []
	ids = {}
	for seq in read_fasta(fastaFile):
		seqInfo = seqPattern.match(seq.description)
		if not seqInfo:
			raise RuntimeError("Sequnces must contain HSP-Tiler headers.")

		contig = seq.name #.split("_")[0]
		gi = seqInfo.group(1)
		score = float(seqInfo.group(2))
		yield contig, gi, score 
		#ids[seq.name.split("_")[0]] = gi

		#originalScores.append()

	#return izip(ids.keys(), ids.values(), originalScores) 

def run_blastx(fastaFile, blastdb, num_threads=3, blastpath=""):
	"""Calculate updated HSP-Tiler scores. The sequences that were used to 
	create the tile are used to filter the blast database in order to remove 
	hits that were filtered during HSP-Tiler.

	Input:
	fastaFile - path to HSP-Tiler output FASTA
	blastdb - path to the BLAST db
	makegilist - bool. Create a new databse filtered by sequenced used to start the tile
	blastx - path to BLASTX executable
	"""

	print fastaFile, blastdb
	
	print "Running BLAST"
	blastFile = "{}.blast.txt".format(fastaFile)

	if not hasBLAST:
		if blastpath and not blastpath.endswith("/"):
			blastpath += "/"
		cmd = "{}blastx -query {} -db {} -max_target_seqs 10 -outfmt 6 -out {} -num_threads {}".format(blastpath, 
			                                                                                           fastaFile, 
			                                                                           				   blastdb,
			                                                                           				   blastFile,
			                                                                           				   num_threads)
		process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
		out, err = process.communicate()

	else:
		#Use Biopython to call BLASTX, faster? safer?
		blastx_cline = NcbiblastxCommandline(query=fastaFile, db=blastdb, max_target_seqs=10,
	                                     	 outfmt=6, num_threads=4, out=blastFile)
		print >> sys.stderr, "Running BLASTX"
		stdout, stderr = blastx_cline()

	return open(blastFile)

def compare(results, original_scores):
	for originalContig, originalGI, originalScore in original_scores:
		try:
			score = max(s for (contig, gi), s in results if contig==originalContig)
		except:
			print "No score for {}, with original score of {}".format(originalContig, originalScore)
			score = originalScore
		yield originalContig, originalScore, score

	"""
	updatedScores = resultsFile
	updatedContig, updatedGI, updatedScore = next(updatedScores)
	originalContig, originalGI, originalScore = next(original_scores)
	done = False
	while not done:
		print "originalContig={}, originalGI={}, originalScore={}".format(originalContig, originalGI, originalScore)
		print "updatedContig={}, updatedGI={}, updatedScore={}".format(updatedContig, updatedGI, updatedScore)
		if originalScore == 0.0:
			#Was not found in the database before or after
			yield originalContig, 0.0, 0.0
			try:
				originalContig, originalGI, originalScore = next(original_scores)
			except StopIteration:
				done = True
			print "next original"
		else:
			if originalContig == updatedContig:
				if originalGI == updatedGI:
					yield updatedContig, originalScore, updatedScore

					while originalContig == updatedContig:
						try:
							updatedContig, updatedGI, updatedScore = next(updatedScores)
						except StopIteration:
							done = True
						print "next updated", updatedContig

					try:
						originalContig, originalGI, originalScore = next(original_scores)
					except StopIteration:
						done = True
					print "next original"

					
				else:
					try:
						updatedContig, updatedGI, updatedScore = next(updatedScores)
					except StopIteration:
						done = True
					print "next updated", updatedContig
			else:
				try:
					updatedContig, updatedGI, updatedScore = next(updatedScores)
				except StopIteration:
					done = True
				print "next updated", updatedContig"""

def read_blastx_bitscores(resultsFile):
	"""Process blastx tab delimeted output. Only saves best hit for each query
	and only save queries who are in the the list of ids

	Input:
	-resultsFile, file-like object of blastx output
	-ids, list of ids of query ID's

	Output:
	-scores, list of bitsocres for each hit
	"""
	results = {}
	for field in csv.reader(resultsFile, delimiter="\t"):
		contig = field[0] #.split("_")[0]
		gi = field[1].split("|")[1]
		updatedScore = float(field[11])
		results[(contig, gi)] = updatedScore
	return results

def write_scores(updated_scores, outfile):
	"""Write scores to File

	Input:
	fasta - file-like object containing FASTA sequences
	original_scores - a list of scores
	updated_scores - a list of scores
	"""
	print >> outfile, "Contig\tOld Score\tNew Score\thsp_tiler_score"
	#for sequence, oldScore, newScore in izip(read_fasta(fasta), original_scores, updated_scores):
	for contig, originalScore, updatedScore in updated_scores:
		print "Writing", contig, originalScore, updatedScore 
		if originalScore != 0.0:
			score = originalScore/updatedScore
		else:
			score = 0.0
		print >> outfile, "{}\t{}\t{}\t{}".format(contig, originalScore, updatedScore, score)

def read_scores(score_file):
	"""Read in a score file written by the write_scores function. Useful for re running script
	to edit graphs or view the in iPython notebook

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

def parse_args():
	parser = argparse.ArgumentParser(description="Score HSP-Tiler output FASTA files")
	parser.add_argument("-f", "--fasta", 
	                    required=True,
	                    type=argparse.FileType('r'),
	                    help="Blastx xml file containing annotations for sequences")
	runOrRead = parser.add_mutually_exclusive_group(required=True)
	runOrRead.add_argument("-b", "--blastdb",
	                       help="Location of blast db to calcuate evalue changes")
	runOrRead.add_argument("-r", "--results",
		                   type=argparse.FileType('r'),
		                   help="Load in a BLAST output file, if BLAST was already run")
	parser.add_argument("--blastpath",
	                    required=False,
	                    default="",
	                    help="Location of blast executibles, if not in user's PATH")

	#Define output
	parser.add_argument("-o", "--outfile",
	                    type=argparse.FileType('wt'),
	                    default=sys.stdout,
	                    help="File to save corrected sequences")
	return parser.parse_args()

if __name__ == "__main__":
	#Parse arguments
	args = parse_args()

	fasta1, fasta2 = tee(args.fasta)

	#Score
	original_scores = get_original_scores(fasta1)

	if args.blastdb is not None:
		blastFile = run_blastx(args.fasta.name, args.blastdb, blastpath=args.blastpath)
		updated_scores = read_blastx_bitscores(blastFile)
	elif args.results is not None:
		print "Here", args.results
		updated_scores = read_blastx_bitscores(args.results)
	else:
		raise RuntimeError("Must run BLASTX (--blastdb) or read in a BLASTX file (--results)")	

	#Compare scores
	scores = compare(updated_scores, original_scores)
	#Save scores
	write_scores(scores, args.outfile)

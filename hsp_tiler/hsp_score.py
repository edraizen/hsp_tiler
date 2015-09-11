#!/usr/bin/python2.7
#Author: Eli Draizen 
#File: hsp_score.py

#Standard Libraries
import sys
import re
import argparse
import csv
import subprocess
from itertools import izip, tee, groupby
from collections import defaultdict

#Required Libraries
from scipy.stats import t, ttest_ind

#Custom imports
from BlastXMLParser import BlastXMLParser
from BlastTabParser import BlastTabParser
from izip_missing import izip_missing
from read_fasta import read_fasta

"""Score the output from hsp_tiler"""

class HSPScore(object):
    """
    """
    def __init__(self, original, updated, normalize=None, trinity=False):
        """
        Parameters:
        ___________
        original : generator that yields Matches objects
        updated : generator that yields Matches objects
        """
        self.original = original
        self.updated = updated
        self.normalize = normalize
        self.original_bitscores = []
        self.updated_bitscores = []
        self.trinity = trinity

        if self.normalize:
            self.max_bitscore = {}
            for matches in self.normalize:
                for match in matches.matches:
                    print matches.contig_name, match.hitID, matches.contig_name.split("|")[1], match.hitID
                    #if (("|" in matches.contig_name and not matches.contig_name.split("|")[1] == match.hitID) or
                    #    (not matches.contig_name == match.hitID)
                    #    ):
                    #    continue
                    if matches.contig_name.split("|")[1] == match.hitID:
                        self.max_bitscore[match.hitID] = match.bitscore
            #self.normalize = {match.hitID: match.bitscore for matches in \
            #    self.normalize for match in matches.matches if \
            #    matches.contig_name.split("|")[1] == match.hitID} #Only works for NCBI formatted headers
            print self.max_bitscore

        else:
            self.updatedIDs = []

    def compare(self):
        """Compare the original scores to the updated scores. The sequence that
        started tile's score is looked at before and after.

        Parameters:
        ___________
        updated_scores : dictionary with key a tuple of contig and gi and value the score, 
            from read_blastx_bitscores
        original_scores : iterable that returns the original contig name, the
            original GI that starte the tile, and the score of the original
            BLAST run.
        """
        def parseTrinity(matches):
            """Parse trintiy headers. Returns KeyA and keyB"""
            parts = matches.contig_name.split("_")
            key = (int(parts[0][4:]), int(parts[1][1:]), int(parts[2][3:]))
            return key, key

        def parseNCBI(matches):
            """Parse NCBI headers. Returns KeyA and keyB"""
            return matches.contig_name, matches.contig_name.split(" ")[0]

        keyA, keyB = parseTrinity if self.trinity else parseNCBI
        for original, updated in izip_missing(self.original, 
                                              self.updated,
                                              #key=lambda x: ,
                                              keyA=keyA, 
                                              keyB=keyB,
                                              fillvalue=None,
                                              verbose=True,
                                              ):

            #print "ORIG, UPDATED", original, updated
            if original is None:
                print "No original"
                continue

            """if original.contig_name == "comp12449_c0_seq4" and updated.contig_name == "comp12449_c0_seq4":
                parse = True

            if not parse:
                print original.contig_name
                continue"""

            
            original_contig = original.contig_name
            original_matches = original.matches
            o_score_normalized = None

            if len(original_matches) == 0:
                print "No matches for {}".format(original_contig)
                continue

            original_match = original_matches[0] #Get the first hit that started the tile

            if self.normalize:
                if original_match.hitID in self.max_bitscore:
                    print "here", self.max_bitscore[original_match.hitID]
                    o_score_normalized = original_match.normalized_bitscore(self.max_bitscore[original_match.hitID])
                else:
                    print "Huh?"
                    #Skip if not in nomalized
                    continue
            else:
                self.original_bitscores.append(original_match.bitscore)
            
            if updated is None:
                print "No score for {}, with original score of {}".format(original_contig, original_match.bitscore)
                updated_score = original_match.bitscore
            else:
                print "WHAt"
                matches = list(updated.matches)
                print matches
                updated_contig = updated.contig_name
                updated_score = 0.0
                u_score_normalized = None
                for match in matches:
                    print match.hitID, original_match.hitID
                    if match.hitID == original_match.hitID and (updated_score == 0.0 or not u_score_normalized):
                        print "Yes score"
                        if self.normalize:
                            u_score_normalized = match.normalized_bitscore(self.max_bitscore[match.hitID])
                        else:
                            updated_score = match.bitscore

                        if self.normalize:
                            #normalized_bitscore = match.normalized_bitscore(self.max_bitscore[match.hitID])
                            
                            if o_score_normalized:
                                print "Adding original normalized!", match.hitID
                                self.updated_bitscores.append(u_score_normalized)
                                self.original_bitscores.append(o_score_normalized)
                                break
                        else:
                            self.updated_bitscores.append(match.bitscore)
                            self.updatedIDs.append(match.hitID)
                            print "Saving {} for {}, from {}".format(match.bitscore, original_contig, original_match.bitscore)
                            break
                if not self.normalize:
                    print original_contig, o_score_normalized, u_score_normalized

                if not self.normalize and updated_score == 0.0:
                    print "No Score for", original_match.hitID
                    updated_score = original_match.bitscore

            if self.normalize:
                if o_score_normalized and u_score_normalized:
                    yield original_contig, o_score_normalized, u_score_normalized,match.hitID
                else:
                    print "Error for", original_contig
            else:
                #print original_contig
                yield original_contig, original_match.bitscore, updated_score, match.hitID

    def mean_bitscore(self, scores):
        """
        """
        print len(scores)
        return sum(score for score in scores)/float(len(scores))

    def write(self, outfile):
        """Write scores to File

        Parameters:
        ___________
        outfile: file-like object to write updated scores
        """
        print >> outfile, "#HSP-Tiler v2.0"
        print >> outfile, "#Contig\tOld Score\tNew Score\thsp_tiler_score"
        #for sequence, oldScore, newScore in izip(read_fasta(fasta), original_scores, updated_scores):
        for contig, originalScore, updatedScore, hitID in self.compare():
            #print "Writing", contig, originalScore, updatedScore 
            if originalScore != 0.0:
                score = updatedScore/originalScore if originalScore else 0.0
            else:
                score = 0.0
            print >> outfile, "{}\t{}\t{}\t{}\t{}".format(contig, originalScore, updatedScore, score, hitID)

        if self.normalize:
            o_bits = self.mean_bitscore(self.original_bitscores)
            u_bits = self.mean_bitscore(self.updated_bitscores)
        else:
            o_bits = self.original_bitscores
            u_bits = self.updated_bitscores
        print o_bits
        print u_bits
        print >> outfile, "#Mean Bitscore (Original):", o_bits
        print >> outfile, "#Mean Bitscore (Updated):", u_bits
        t_test = ttest_ind(self.original_bitscores, self.updated_bitscores)
        print >> outfile, "#Two sample t-test (Same Variance): t =", t_test[1], "p ="
        print >> outfile, "#Satisifies 95% C.I.:", t_test[1] <= 0.05
        print >> outfile, "#Satisifies 99% C.I.", t_test[1] <= 0.01
        t_test = ttest_ind(self.original_bitscores, self.updated_bitscores, equal_var=False)
        print >> outfile, "#Two sample t-test (Diff Variance): t =", t_test[1], "p ="
        print >> outfile, "#Satisifies 95% C.I.:", t_test[1] <= 0.05
        print >> outfile, "#Satisifies 99% C.I.", t_test[1] <= 0.01

def read_scores(score_file):
    """Read in a score file written by the write_scores function. Useful for re running script
    to edit graphs or view the in iPython notebook

    Parameters:
    score_file, a file-like object containing tab delimeted score informations

    Output:
    Three lists, old scores, new scores, and the hsp-tile, score
    """
    def no_comments(it):
        for line in it:
            if not line.startswith("#"):
                yield line

    old_scores = []
    new_scores = []
    scores = []
    for fields in csv.reader(no_comments(score_file), delimiter="\t"):
        old = float(fields[1])
        new = float(fields[2])
        score = float(fields[3])
        yield old, new, score

def parse_args():
    parser = argparse.ArgumentParser(description="Score HSP-Tiler output FASTA files")
    parser.add_argument("--original", 
                        required=True,
                        type=argparse.FileType('r'),
                        help="Blastx output file containing annotations for original sequences")
    parser.add_argument("-r", "--results",
                        type=argparse.FileType('r'),
                        help="Blastx output file containing annotations for updated sequences")
    parser.add_argument("-n", "--normalize",
                        default=None,
                        type=argparse.FileType('r'),
                        help="Blast output file containing annotations for each hit from the updated sequences against themselves")
    parser.add_argument("--oformat", 
                        default="blastxml", 
                        choices=("blastxml", "blasttab", "pslx", "rs2"),
                        help="Format of the original annotation file. Can use output from BLASTX (xml (default), or tabs with specified format), BLAT PSLX, and RapSearch2.")
    parser.add_argument("--uformat", 
                        default="blasttab", 
                        choices=("blastxml", "blasttab", "pslx", "rs2"),
                        help="Format of the updated annotation file. Can use output from BLASTX (xml (default), or tabs with specified format), BLAT PSLX, and RapSearch2.")
    parser.add_argument("--nformat", 
                        default="blasttab", 
                        choices=("blastxml", "blasttab", "pslx", "rs2"),
                        help="Format of the normaliziation annotation file. Can use output from BLASTX (xml (default), or tabs with specified format), BLAT PSLX, and RapSearch2.")
    parser.add_argument("--trinity",
                        required=False, 
                        default=False, 
                        action="store_true")
    #evlaue cutoff
    parser.add_argument("-e", "--evalue_cutoff",
                        type=float,
                        default=1e-10,
                        help="Only allow Blast Hist less than or equal to cutoff. Default is 1e-10")
    #Use every hit for eqch query instead of just the 1st one
    parser.add_argument("--useHit1",
                        default=False,
                        action="store_true",
                        help="Use only the first, high scoring hit for each sequence. Defualt is to use all hits from annotation. Default is to all hits. Optional.")
    #family/species filter
    parser.add_argument("--filter",
                        required=False,
                        help="Filter out blastx hits from a given species and fmily name using full scientific names only for now. Must be the same as the one used in HSP-Tiler")
    parser.add_argument("--filterType",
                        default=0, #0 for ncbi, 1 for uniprot
                        help="Regex to retreive taxon info from hit. Use 0 for NCBI, or 1 for uniprot, or the full regex query")

    #Define output
    parser.add_argument("-o", "--outfile",
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="File to save corrected sequences")
    return parser.parse_args()

if __name__ == "__main__":
    #Parse arguments
    args = parse_args()

    #Score
    if args.oformat == "blastxml":
        original_scores = BlastXMLParser(args.original, 
                                         allHits=not args.useHit1, 
                                         evalue=args.evalue_cutoff, 
                                         taxFilter=args.filter, 
                                         taxFilterType=args.filterType)
    elif args.oformat == "blasttab":
        original_scores = BlastTabParser(args.original, 
                                         allHits=not args.useHit1, 
                                         evalue=args.evalue_cutoff, 
                                         taxFilter=args.filter, 
                                         taxFilterType=args.filterType,
                                         qframe=False)
    else:
        raise RuntimeError("Invalid original format")

    if args.uformat == "blastxml":
        updated_scores = BlastXMLParser(args.results)
        #print list(updated_scores.parse().next().contig_name)
    elif args.uformat == "blasttab":
        updated_scores = BlastTabParser(args.results, qframe=False, query_seq=False)
    else:
        raise RuntimeError("Invalid updated format")

    if args.normalize:
        normalized_scores = BlastTabParser(args.normalize, qframe=False, query_seq=False).parse()
    else:
        normalized_scores = None

    #Compare scores
    scores = HSPScore(original_scores.parse(), updated_scores.parse(), normalize=normalized_scores, trinity=args.trinity)
    #Save scores
    scores.write(args.outfile)

    if not args.normalize:
        with open("hitIDs.txt", "w") as uID:
            for i in scores.updatedIDs:
                print >> uID, i


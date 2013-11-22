#!/usr/local/bin/python

#######################################
#Imports
#######################################

#Standard Libraries
import argparse
import os, sys
import re
import string
from collections import deque
from itertools import izip, tee
import subprocess

#Custom Libraries
from read_fasta import read_fasta
from BlastXMLParser import BlastXMLParser

#######################################
#Classes
#######################################

#HSPs
class Hsp(object):
    """Holds information about a given High-scoring
    Sequence Pair (HSP) returned from BlastX. The 
    information held is the strand of the query 
    that the HSP is on, the frame it is in relative 
    to the query sequence, the start and end positions 
    of the HSP relative to the query sequence and the 
    start and end positions relative to the protein 
    sequence of the hit, the e-value and bitscore of 
    the hit, the amino acid sequence of the HSP and the 
    number of the HSP. This object also keep track to see
    if HSP has already be tiled.
    """
    def __init__(self, hsp, contig):
        """Initialise an Hsp object, with the XML 
        representation of the hsp, and a contig.  

        Input:
        hsp - Element object containing HSP information
        contig - Seqeunce object
        """
        # 1 = frames 1,2,3 -1 = frames -1,-2,-3 
        if 3 >= int(hsp.find('Hsp_query-frame').text) > 0:
            self.strand = 1
        elif 0 > int(hsp.find('Hsp_query-frame').text) >= -3:
            self.strand = -1
        else: # Something is wrong - exit here!
            raise SystemExit ('Cannot parse HSP data correctly')
        self.frame = abs(int(hsp.find('Hsp_query-frame').text))
        self.query_start = int(hsp.find('Hsp_query-from').text) #position in nts relative to query
        self.query_end = int(hsp.find('Hsp_query-to').text)
        self.hit_start = int(hsp.find('Hsp_hit-from').text) #position in aas relative to hit
        self.hit_end = int(hsp.find('Hsp_hit-to').text)
        self.evalue = float(hsp.find('Hsp_evalue').text)
        self.score = float(hsp.find('Hsp_score').text)
        self.aa_seq = self.remove_x(hsp,contig) 
        self.used = False # True marks hsps that have already been incorporated into the tile
        self.num = int(hsp.find('Hsp_num').text)

    def remove_x (self, hsp, contig):
        """Remove Xs and gaps from hsp aa sequence by 
        substituting from contig translation to fix sequences 
        that have repeptative or low-complexity regions, or 
        filtering was turned on in BlastX. This method creates 
        regular expression from the HSP protein sequence. 
        If the HSP protein sequence were FGTPPXPYII, the 
        regular expression created would be FGT....YII. This 
        is used to search all the translations of the contig 
        or a matching sequence.

        Input:
        hsp - Element object containing HSP information
        contig - Seqeunce object
        """
        no_x = re.sub(r'-', r'', hsp.find('Hsp_qseq').text) # remove gaps 
        start = re.sub(r'[X|x]', r'.', no_x[:3]) # first 3 aas - replace X with . for re matching
        start = re.sub(r'\*', r'.', start) # change * to _ to match stop codons correctly
        end = re.sub(r'[X|x]', r'.', no_x[-3:]) # last 3 aas
        end = re.sub(r'\*', r'.', end)
        matches = self.match(len(no_x), contig, start, end)
        if len(matches) != 1: # if there are no matches, try looking for unknown bases as well
            start= "({}|J)({}|J)({}|J)".format(*start)
            end = "({}|J)({}|J)({}|J)".format(*end)
            matches = self.match(no_x, contig, start, end)
        match = matches[0].group() # and assuming there is only one match, the index is correct
        if len(match) != len(no_x): # sanity check
            raise SystemExit("Cannot match HSP {} to Contig {}".format(hsp.num, contig.name))
        else:
            while 'X' in no_x:
                xs = re.search(r'X+', no_x) # search for X or runs of X
                s = xs.start() # indices of xs matched above
                e = xs.end()
                subst = match[s:e] # portion of contig to use to replace Xs
                no_x = re.sub(xs.group(), subst, no_x, count=1) # only replace first instance in each loop
        no_x = re.sub(r'J', r'X', no_x) # replace J (from unknown translations) with X as these are genuine unknowns
        return no_x

    def match(self, len_seq, contig, start, end):
        """Find the protein sequence encoded by the contig 
        that matches the hsp query sequence. Uses the first
        three and last three aas, separted by any character.
        from the query to build a regular extansion,

        Input:
        len_seq - int, Length of sequence that has no Xs or gaps
        contig - Sequence object of contig
        start - First 3 aas of query sequence
        end - Last 3 aas of query sequence 
        """
        anychar_pattern = "."*(len_seq-6)
        pattern = r"{}{}{}".format(start, anychar_pattern, end) # regex for matching
        matches = []
        for sequence in contig.sixframe.keys():
            match = re.search(pattern, contig.sixframe[sequence]) # search all protein translations for regex
            if match != None:
                matches.append(match)
        return matches

    def printer (self):
        """Writes hsp information to global log File
        """
        print >> logfile, "Hsp {}".format(self.num)
        print >> logfile, self.aa_seq
        print >> logfile, "Strand: {}".format(self.strand)
        print >> logfile, "Frame: {}".format(self.frame)
        print >> logfile, "Start: {}, End: {}".format(self.query_start, self.query_end)

#######################################

class Tile_Path(object):
    """An object to hold information about the tile and 
    methods for adding new HSPs to the tile. If an HSP 
    is added 5’ of the tile, it is appended to the left 
    of the queue and if it is added 3’ of the tile, it is 
    appended to the right of the list. In this way, the 
    order of the HSPs in the list represents the order of 
    HSPs in the tile form 5’ to 3’. The nucleotide 
    sequence of the tile is derived by slicing the 
    appropriate portion of the contig sequence using the 
    indices for the start and end positions of the HSPs.
    """
    def __init__(self, best_hsp, contig): 
        """Initialise the tile with a single HSP to signify the start 
        and end positions of the tile relative to the contig and the 
        protein sequence of the tile.

        Input:
        best_hsp - Hsp object of the highest scoring Hsp
        contig - Sequence object of the parent contig
        """
        #Positions are of tiled construct relative to contig (query)
        self.contig = contig
        self.start = best_hsp.query_start
        self.end = best_hsp.query_end
        self.aa_seq = best_hsp.aa_seq 
        self.nt_seq = self.contig.sequence[self.start-1 : self.end] # nt sequence of tile from position on contig
        self.hsps = deque([best_hsp]) # double-ended queue allows insertion of hsps at either end
        self.strand = best_hsp.strand
        self.frame = best_hsp.frame
        self.added = True # Has anything been added this iteration? Initialise as true to allow loop to start! 
        self.tile = False # Has anything other than the first hsp been added to the tile?

    def __str__(self):
        """Print tile sequence if changes have been made, 
        else return unchanged contig.
        """
        if self.tile:
            seq = self.nt_seq
        else:
            seq = self.contig.sequence

        tmp = ">{}\n".format(self.contig.name)
        for i in range(0, len(seq), 60):
            tmp += "{}\n".format(seq[i:i+60])
        return tmp[:-1]

    def add_hsp(self, hsp, location, fill=0): 
        """Generic method to add an hsp to the tile path. 
        This method adjusts the indices of the HSP using 
        the value of the fill and finds the nucleotide 
        sequence created by slicing the contig sequence 
        using those indices. This nucleotide sequence is 
        concatenated with the pre-existing nucleotide 
        sequence of the tile and the new amino acid 
        sequence created by translating the nucleotide 
        sequence. The incoming HSP is also added to the 
        appropriate end of the double-ended queue.

        Input:
        hsp - Hsp object if Hsp to add.
        location - either "upstream" or "downstream." The 
            location of the HSP relative to the existing tile.
        fill - int. the absolute difference between the end 
            of the tile and the end of the HSP; adjust to mend
            possible frameshifts
        """
        if location == 'upstream': # add hsp to 5' end
            self.start = hsp.query_start # adjust nucleotide coordinates relative to contig
            self.hsps.appendleft(hsp) # append hsp to list of those included at correct end
            self.nt_seq = self.contig.sequence[self.start-1: hsp.query_end+fill] + self.nt_seq
        elif location == 'downstream': # add hsp to 3' end
            self.end = hsp.query_end
            self.hsps.append(hsp)
            self.nt_seq = self.nt_seq + self.contig.sequence[hsp.query_start-fill-1: self.end]
        else: # something is wrong, exit here
            raise SystemExit ('Error: Cannot add hsp to tile')
        self.aa_seq = translate_sequence(self.nt_seq, self.strand)
        hsp.used = True # mark hsp as used
        self.added = True # something has been incorporated into tile
        self.tile = True

    def add_overlap(self, hsp, location, distance):
        """Add an hsp where it overlaps the existing tile.
        Where the HSP overlaps the contig, the fill is 
        subtracted from rather than added to the incoming HSP.
        This trims nucleotides from the sequence being added to 
        the tile until they translate in the same frame.

        Input:
        hsp - Hsp object if Hsp to add.
        location - either "upstream" or "downstream." The 
            location of the HSP relative to the existing tile.
        distance - distance between hsp and this tile, used as fill value
        """
        print >> logfile, "Hsp number {} is {} of the tile and overlapping by {} nucleotides".format(hsp.num, location, distance+1)
        if (distance +1) % 3 == 0:
            self.add_hsp(hsp, location, 0-(distance+1))
        else:
            print >> logfile "Hsp number {} not in same frame as tile - correcting".format(hsp.num)
            while (distance +1) % 3 !=0:
                distance -= 1
            self.add_hsp(hsp, location, 0-(distance+1))
        print >> logfile, "Overlapping HSP number {} added {} of tile".format(hsp.num, location)
        self.printer()

    def add_insert(self, hsp, location, distance):
        """Add an hsp where there is an insertion in the contig 
        between the current hsp and the existing tile. When the 
        incoming HSP is not in the same frame as the tile, it is 
        corrected. 

        Input:
        hsp - Hsp object if Hsp to add.
        location - either "upstream" or "downstream." The 
            location of the HSP relative to the existing tile.
        distance - distance between hsp and this tile, used as fill value
        """
        difference = distance -1
        if difference == 1 or difference == 2: # assume incorrect insertion - delete inserted nts
            self.add_hsp(hsp, location)
        elif difference % 3 == 0: # assume correct insertion - include nts and translate
            self.add_hsp(hsp, location, difference)
        else: # insertion with frameshift - attempt to fix
            self.correct_frame(hsp, location, distance) 
        
        print >> logfile, "Insertion in contig relative to hit"
        print >> logfile, "Hsp number {} added {} of tile".format(hsp.num, location)
        self.printer()

    def add_subst(self, hsp, location, distance):
        """Add an hsp where there is a substitution or deletion 
        in the contig between the current hsp and the existing tile.
        When the incoming HSP is not in the same frame as the tile, 
        it is corrected. 

        Input:
        hsp - Hsp object if Hsp to add.
        location - either "upstream" or "downstream." The 
            location of the HSP relative to the existing tile.
        distance - distance between hsp and this tile, used as fill value
        """
        difference = distance -1
        if difference == 0 or difference % 3 == 0: # deletion or substitution respectively - include and translate nts if present
            self.add_hsp(hsp, location, difference)
        else: # substitution with frameshift
            self.correct_frame(hsp, location, distance)

        print >> logfile, "Deletion or substitution in contig relative to hit"
        print >> logfile, "Hsp number {} added {} of tile".format(hsp.num, location)
        self.printer()

    def correct_frame(self, hsp, location, distance):
        """Make a correction between hsps when they are not in the same frame.
        Frames are corrected by subtracting from the difference 
        until it is a multiple of three before adding it to the
        tile. This trims nucleotides from the sequence being 
        added to the tile until they translate in the same frame.

        Input:
        hsp - Hsp object if Hsp to add.
        location - either "upstream" or "downstream." The 
            location of the HSP relative to the existing tile.
        distance - distance between hsp and this tile, used as fill value
        """
        difference = distance -1
        while difference % 3 != 0: # remove nts until difference is divisible by 3
            difference = difference -1
        self.add_hsp(hsp, location, difference) # add hsp with the cropped nts and translate

        print >> logfile, "Hsp number {} not in same frame as tile - correcting".format(hsp.num)

    def printer(self):
        """Writes hsp information to global log File
        """
        print >> logfile, "Sequence: {}".format(self.nt_seq)
        print >> logfile, "Start: {}, End: {}".format(self.start, self.end)
        print >> logfile, "Protein: {}".format(self.aa_seq)

#######################################
#Global functions
#######################################

def translate_sequence(sequence, strand=1):
    """Translates and sequence in the first frame, 
    forward (default) or reverse.

    Input:
    sequence - String contianing sequence to translate
    strand - 1 for forward (default); -1 for reverse
    """
    if strand == -1:
        sequence = revcomp(sequence)
    protein = codons (sequence)
    return protein

gencode = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
def codons (sequence):
    """Return a protein string for a sequence in frame 1

    Input:
    sequence - String contianing sequence to translate
    """
    aas = []
    codons = [sequence[i:i+3] for i in range (0, len(sequence), 3)]
    for codon in codons:
        if len(codon) == 3: # 'codons' that aren't 3 bases should only occur at the end of a sequence - can be discarded
            if codon.upper() in gencode:
                aas.append(gencode[codon.upper()])
            else: 
                aas.append('J') # substitute J if triplet can't be translated - X not used to prevent unending loops in replace_x
    protein = ''.join(aas)
    return protein

complement_table = string.maketrans("ACGT", "TGCA")
def revcomp(dna):
    """returns the reverse complement of a sequence

    Input:
    dna - string containg dna sequence
    """
    return dna[::-1].translate(complement_table) 

def sixframe(dna):
    rev_comp = revcomp(dna)
    sixframe = {}
    for i in range (0, 3):
        seq = dna[i:] # get frames positive strand
        protein = codons(seq)
        sixframe [i+1] = protein
        rev = rev_comp[i:] # reverse complement and do again to get the other 3 frames
        protein = codons(rev)
        sixframe [0-i-1] = protein
    return sixframe

def parse_args():
	# Parsing command line options
	parser = argparse.ArgumentParser(description="Takes a fasta file of sequences and a BLASTX annotation of that file in xml format.  Attempts to tile Hsps for the highest scoring hit for each sequence, correcting frameshifts in order to improve subsequent annotations.")
	# name of fasta file 
	parser.add_argument("-f", "--fasta", 
	                    required=True, 
	                    type=argparse.FileType('r'),
	                    help="Fasta file containing sequences")
	# name of annotation file (in xml format as code currently stands)
	parser.add_argument("-a", "--annotation", 
	                    required=True, 
	                    type=argparse.FileType('r'),
	                    help="Blastx xml file containing annotations for sequences")
	# gap limit
	parser.add_argument("-g", "--gap_limit", 
	                    type=int, 
	                    default=15, 
	                    help="Cutoff for distance between hsps. If the gap between hsps in nucleotides is greater than this value, the hsp will not be added to the tile.  Default = 15nt")
	#family/species filter
	parser.add_argument("--allHits",
						default=False,
						help="Use all hits from annotation. Default is to use only the first, high scoring hit. Optional.",
						action="store_true")
	parser.add_argument("--filter",
	                    required=False,
	                    help="Filter out blastx hits from a given species and fmily name using full scientific names only for now")
	#Define output
	parser.add_argument("-o", "--outfile",
	                    type=argparse.FileType('wt'),
	                    default="hsp_tiler_output.fa",
	                    help="File to save corrected sequences")
	parser.add_argument("-l", "--logfile",
	                    type=argparse.FileType('wt'),
	                    default=sys.stdout,
	                    help="File to save log")

	# print help message if no arguments are given
	if len (sys.argv) == 1:
		print "here?" 
		parser.print_help()
		sys.exit(1)

	#Parse args
	return parser.parse_args()

def run(fasta, annotation, gap_limit, allHits, filter, outfile):
	#Initilise BLAST parser
	parser = BlastXMLParser(annotation, allHits=allHits, taxFilter=filter)

	# get list of hsps for highest scoring hit for each contig - go to next contig if no hits
	for contig, queryID in izip(read_fasta(fasta), parser.parseQuery()):
	    if not contig.name.strip() == queryID:
	        raise RuntimeError("BLAST Queries are not in the same order as FASTA file: {}, {}".format(contig.name.strip(), 
	        																						  queryID))

	    print >> logfile, "\nNext contig:"
	    print >> logfile, str(contig)

	    contig.sixframe = sixframe(contig.sequence)

	    
	    #Add each hsp to list
	    hsp_list = [Hsp(hsp, contig) in parser.parseHsp()]
	    
	    if not hsp_list:
	        print >> logfile, "No hits for {}".format(contig.name)
	        print >> outfile, str(contig)
	        continue
	         
	    # Initialise tile using highest scoring hsp (first in list) 
	    tile = Tile_Path (hsp_list[0], contig) 
	    hsp_list[0].used = True # mark that this hsp has been used
	    print >> logfile, "Starting tile..."
	    hsp_list[0].printer()
	    tile.printer()

	    # go to next contig if only one hsp in list
	    if len(hsp_list) == 1:
	        print >> logfile, "Only one hsp for this hit.  Tile is complete"
	        print >> outfile, str(contig)
	        continue # next contig

	    # if >1 hsp, attempt to create a tile
	    else: 
	        # run loop until no hsps were added in the previous iteration
	        while tile.added == True: 
	            tile.added = False # loop will exit if this is not set to True by calling add.hsp() in this iteration

	            # iterate hsps for this contig
	            for hsp in hsp_list:

	                # ignore hsps that are already in the tile
	                if hsp.used == False: 
	                    logfile.write('\nNext hsp...\n')
	                    hsp.printer()

	                    # ignore hsps that are not on the same strand as the tile
	                    if hsp.strand != tile.strand: 
	                        print >> logfile, "Hsp number {} is not on the same strand as the tile - cannot tile".format(hsp.num)
	                        hsp.used = True
	                        continue # next hsp

	                    # ignore hsps that are completely within the tile
	                    if tile.start <= hsp.query_start and hsp.query_end <= tile.end: 
	                        print >> logfile, "Hsp {} completely within tile -skipped".format(hsp.num)
	                        hsp.used = True # mark as used but do not add to tile
	                        tile.printer()
	                        continue

	                    # determine whether incoming hsp is 5' or 3' of tile
	                    if hsp.query_start <= tile.start: # hsp nt seq 5' of tile
	                        location = 'upstream'
	                        distance = abs(hsp.query_end - tile.start)
	                    elif hsp.query_end >= tile.end: # hsp nt seq 3' of tile
	                        location = 'downstream'
	                        distance = abs(hsp.query_start - tile.end)

	                    # determine whether hsps overlap and add if they do
	                    if hsp.query_start <= tile.start <= hsp.query_end or hsp.query_start <= tile.end <= hsp.query_end:
	                        tile.add_overlap(hsp, location, distance) 

	                    # otherwise they do not overlap
	                    else: 

	                        # check hsp locations make sense before continuing
	                        if location == 'upstream': 
	                            assert hsp.query_end < tile.start 
	                        elif location == 'downstream':
	                            assert hsp.query_start > tile.end 

	                        # check if hsp is close enough to tile to incorporate
	                        if distance > gap_limit: 
	                            print >> logfile, "Too far away to tile"
	                        # go to next hsp but do not mark as incorporated as may be able to add in a later iteration    
	                            continue 

	                        # if contig is close enough to attempt tiling                
	                        else: 
	                            print >> logfile"Hsp {} is {} nucleotides {} of tile".format(hsp.num), str(distance -1), location)
	                            assert distance -1 >=0, 'Cannot tile hsps' #sanity check

								#Check whether the gap between the hsps relative to the HIT is 0.  
								#If this is true but there is a nucleotide gap relative to the contig, 
								#something must have been inserted into the contig.  Statement is unwieldy 
								#because test must be made for every combination of hsps upstream and 
								#downstream of the tile and on positive and negative strands.
	                            if (location == 'upstream' and ((tile.strand == 1 and tile.hsps[0].hit_start - hsp.hit_end -1 == 0) or (tile.strand == -1 and hsp.hit_start - tile.hsps[0].hit_end -1 ==0))) or (location == 'downstream' and ((tile.strand == 1 and hsp.hit_start - tile.hsps[-1].hit_end -1 == 0) or (tile.strand == -1 and tile.hsps[0].hit_start - hsp.hit_end -1 == 0))): 
	                                tile.add_insert(hsp, location, distance)

								#If the above is not true, check whether the gap between the hsps 
								#relative to the hit is > 0.  If this is true and there is no nucleotide 
								#gap relative to the contig, something must have been deleted from the contig.  
								#If this is true and there is a nucleotide gap, the nucleotides in that gap 
								#must have been substituted in the contig. Statement is unwieldy as above
	                            elif (location == 'upstream' and ((tile.strand == 1 and tile.hsps[0].hit_start - hsp.hit_end -1 > 0) or (tile.strand == -1 and hsp.hit_start - tile.hsps[0].hit_end -1 > 0))) or (location == 'downstream' and ((tile.strand == 1 and hsp.hit_start - tile.hsps[-1].hit_end -1 > 0) or (tile.strand == -1 and tile.hsps[0].hit_end - hsp.hit_start -1 > 0))): 
	                                tile.add_subst(hsp, location, distance)
	                           
								#as we have already dealt with hsps that overlap the tile, if neither 
								#of the above are true, there is a problem somewhere
	                            else:
	                                print >> logfile, "Contig {}: HSP{} could not be tiled".format(contig.name, hsp.num)
	                                hsp.used = True
	                                continue

	    # add unchanged contig to output if no changes made.  Add tile sequence if changes have been made.
	    print >> outfile, tile

#######################################
#Main
#######################################

if __name__ == "__main__":
	#Parse arguments
	args = parse_args()

	logfile = args.logfile

	#Run hsp_tiler
	run(args.fasta, 
		args.annotation, 
		args.gap_limit, 
		args.allHits, 
		args.filter, 
		args.outfile)

	#Cleanup
	logfile.close()
	args.outfile.close()

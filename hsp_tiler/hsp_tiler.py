#!/usr/local/bin/python
# Author: Eli Draizen, Kathryn Crouch
# Date: 16-3-2014
# File: hsp_tiler.py

#######################################
#Imports
#######################################

#Standard Libraries
import argparse
import sys
import re
import string
from collections import deque
from itertools import product, tee, imap
from datetime import datetime
from multiprocessing import Pool

#Custom Libraries
from read_fasta import read_fasta
from BlastXMLParser import BlastXMLParser
from BlatPSLXParser import BlatPSLXParser
from Sequence import Sequence, DNASequence
from count_kmers import count_kmers, IUPAC_N
from read_codon_usage import read_codon_frequencies

#######################################
#Global Variables
#######################################

logfile = sys.stderr

UPSTREAM = 0
DOWNSTREAM = 1

#######################################
#Classes
#######################################

class Tile_Path(object):
    """An object to hold information about the tile and methods for adding new 
    HSPs (Matches) to the tile. If an HSP is added 5' of the tile, it is 
    appended to the left of the queue and if it is added 3' of the tile, it is 
    appended to the right of the list. In this way, the order of the HSPs in 
    the list represents the order of HSPs in the tile form 5' to 3'. The 
    nucleotide sequence of the tile is derived by slicing the appropriate 
    portion of the contig sequence using the indices for the start and end 
    positions of the HSPs.

    If best_hsp is None, then there were no hits.
    """
    def __init__(self, best_hsp, contig, saveProteinCoding=True, protein=False, codon_usage=None): 
        """Initialise the tile with a single HSP to signify the start 
        and end positions of the tile relative to the contig and the 
        protein sequence of the tile.

        Parameters:
        ___________
        best_hsp : Mastch object of the highest scoring Match (Hsp)
        contig : Sequence object of the parent contig
        saveProteinCoding : bool. Extend sequence to cover entire protein coding region
        protein : bool. output amino acid sequence
        codon_usage : a dictionary containg the codon usage, must be created 
            using read_codon_usage
        """
        self.contig = contig

        #Positions are of tiled construct relative to contig (query)
        self.start = best_hsp.query_start
        self.end = best_hsp.query_end

        # Get the aa seq of the best Match 
        self.aa_seq = Sequence(name=self.contig.name,
                               sequence=self._getCorrectedAA(best_hsp.query_seq))

        # nt sequence of tile from position on contig
        self.nt_seq = DNASequence(name=self.contig.name,
                                  sequence=self.contig.sequence[self.start:self.end])

        # double-ended queue allows insertion of hsps at either end
        self.hsps = deque([best_hsp]) 

        self.strand = best_hsp.strand
        self.frame = best_hsp.frame
        self.codon_usage = codon_usage
        self.score = best_hsp.score
        self.hitID = best_hsp.hitID
        self.bitscore = best_hsp.bitscore
        self.num = best_hsp.num

        #Has anything been added this iteration? Initialise as true to allow 
        #loop to start!
        self.added = True  

        # Has anything other than the first hsp been added to the tile?
        self.tile = False
        
        #Expand tile to cover entire protein coding region
        self.saveProteinCoding = saveProteinCoding

        #Output proteins sequences instead ot nts
        self.protein = protein

    def __str__(self):
        """Print tile sequence if changes have been made, else return unchanged 
        contig.
        """
        seq = self.nt_seq if not self.protein else self.aa_seq
        seq.description = self._getDescription()

        return str(seq)

    def write(self, outfile):
        """Save tile to file

        Parameters:
        ___________
        outfile : file-like object to write tile

        Return:
        ________
        tile append to outfile
        """
        #Save the protein coding region (tile expanded to closest start and stop codons)
        self.extendReadingFrame()

        if self.protein:
            tile.outputProtein()
        
        #Save the tiles region
        print >> outfile, self

    def outputProtein(self, protein=True):
        """Output protein sequence when being printed

        Paramters:
        __________
        protein : bool. Enable protein output. Default is True.
        """
        self.protein = protein

    def _getDescription(self):
        """Return an updated description with start position of corrected
        contig, frame the corrected contig was in (now in frame 1), the best
        hit that started the tile, and the bit score of that hit. The last two
        options are only available if the annotation file was creating the NCBI
        header.
        """
        return "[Start={};Frame={};GI={};S={}]".format(self.start, 
                                                       self.frame, 
                                                       self.hitID, 
                                                       self.bitscore)

    def add_hsp(self, hsp, location, fill=0): 
        """Generic method to add an hsp to the tile path. This method adjusts 
        the indices of the HSP using the value of the fill and finds the 
        nucleotide sequence created by slicing the contig sequence using those 
        indices. This nucleotide sequence is concatenated with the pre-existing 
        nucleotide sequence of the tile and the new amino acid sequence created 
        by translating the nucleotide sequence. The incoming HSP is also added 
        to the appropriate end of the double-ended queue.

        Parameters:
        ___________
        hsp : Hsp object
        location : either upstream (0) or downstream (1). The location of the 
            Match relative to the existing tile.
        fill : int. the absolute difference between the end of the tile and the 
            end of the HSP; adjust to mend possible frameshifts
        """
        if location == UPSTREAM: 
            # add hsp to 5' end

            # adjust nucleotide coordinates relative to contig
            self.start = hsp.query_start

            # append hsp to list of those included at correct end
            self.hsps.appendleft(hsp)
            
            if hsp.query_end+fill > hsp.query_end:
                #Add gap
                if self.codon_usage:
                    gap = "X"*len(gap)
                else:
                    gap = self.contig.sequence[hsp.query_end:hsp.query_end+fill]
                new_nt = self.contig.sequence[self.start:hsp.query_end]
            else:
                #No Gap
                gap = ""
                new_nt = self.contig.sequence[self.start:hsp.query_end+fill]

            print >> logfile, "Adding HSP upstream"
            print >> logfile, "Original sequence:", self.nt_seq
            print >> logfile, "Fill:", fill, "Gap:", gap
            print >> logfile, "5' Extenison:", new_nt

            self.nt_seq.sequence = "{}{}{}".format(new_nt, gap, self.nt_seq.sequence)

        elif location == DOWNSTREAM: 
            # add hsp to 3' end
            self.end = hsp.query_end
            self.hsps.append(hsp)

            if hsp.query_start-fill < hsp.query_start:
                #Add gap
                if self.codon_usage:
                    gap = "X"*len(fill)
                else:
                    gap = self.contig.sequence[hsp.query_start-fill:hsp.query_start]
                new_nt = self.contig.sequence[hsp.query_start:self.end]
            else:
                #No Gaps
                gap = ""
                new_nt = self.contig.sequence[hsp.query_start-fill:self.end]

            print >> logfile, "Adding HSP downstream"
            print >> logfile, "Original sequence:", self.nt_seq
            print >> logfile, "Fill:", fill, "Gap:", gap
            print >> logfile, "3' Extenison:", new_nt

            self.nt_seq.sequence = "{}{}{}".format(self.nt_seq.sequence, gap, new_nt)

        else: 
            # something is wrong, exit here
            raise RuntimeError('Error: Cannot add hsp to tile')

        self.aa_seq = self.nt_seq.translate_sequence(strand=self.strand)
        
        # mark hsp as used
        hsp.used = True 

        # something has been incorporated into tile
        self.added = True 
        self.tile = True

    def add_overlap(self, hsp, location, distance):
        """Add an hsp where it overlaps the existing tile. Where the HSP 
        overlaps the contig, the fill is subtracted from rather than added to 
        the incoming HSP. This trims nucleotides from the sequence being added 
        to the tile until they translate in the same frame.

        Parameters:
        ___________
        hsp : Hsp object if Hsp to add.
        location : either upstream (0) or downstream (1). The location of the 
            Match relative to the existing tile.
        distance : distance between hsp and this tile, used as fill value
        """
        print >> logfile, "Hsp number {} is {} of the tile and overlapping by",
        print >> logfile,  "{} nucleotides".format(hsp.num, location, distance+1)

        if (distance) % 3 == 0:
            self.add_hsp(hsp, location, 0-(distance))
        else:
            print >> logfile, "Hsp number {} not in same frame as tile - correcting".format(hsp.num)
            while (distance +1) % 3 !=0:
                distance -= 1
            self.add_hsp(hsp, location, 0-(distance))
        print >> logfile, "Overlapping HSP number {} added {} of tile".format(hsp.num, location)
        self.printer()

    def add_insert(self, hsp, location, distance):
        """Add an hsp where there is an insertion in the contig between the 
        current hsp and the existing tile. When the incoming HSP is not in the 
        same frame as the tile, it is corrected. 

        Parameters:
        ___________
        hsp : Hsp object if Hsp to add.
        location : either upstream (0) or downstream (1). The location of the 
            Match relative to the existing tile.
        distance : distance between hsp and this tile, used as fill value
        """

        difference = distance
        if 0 < difference < 3: 
            # assume incorrect insertion - delete inserted nts
            print >> logfile, "Correct Insertion"
            self.add_hsp(hsp, location)
        elif difference % 3 == 0: 
            # assume correct insertion - include nts and translate
            self.add_hsp(hsp, location, difference)
            print >> logfile, "incorrect Insertion"
        else: 
            # insertion with frameshift - attempt to fix
            print >> logfile, "Insertion with frameshift"
            self.correct_frame(hsp, location, distance) 
        
        print >> logfile, "Insertion in contig relative to hit"
        print >> logfile, "Hsp number {} added {} of tile".format(hsp.num, location)
        self.printer()

    def add_subst(self, hsp, location, distance):
        """Add an hsp where there is a substitution or deletion in the contig 
        between the current hsp and the existing tile. When the incoming HSP is
        not in the same frame as the tile, it is corrected. 

        Parameters:
        ___________
        hsp : Hsp object if Hsp to add.
        location : either upstream (0) or downstream (1). The location of the 
            Match relative to the existing tile.
        distance - distance between hsp and this tile, used as fill value
        """
        difference = distance
        if difference == 0 or difference % 3 == 0: 
            # deletion or substitution respectively - include and translate nts if present
            self.add_hsp(hsp, location, difference)
        else: 
            # substitution with frameshift
            self.correct_frame(hsp, location, distance)

        print >> logfile, "Deletion or substitution in contig relative to hit"
        print >> logfile, "Hsp number {} added {} of tile".format(hsp.num, location)
        self.printer()

    def correct_frame(self, hsp, location, distance):
        """Make a correction between hsps when they are not in the same frame.
        Frames are corrected by subtracting from the difference until it is a 
        multiple of three before adding it to the tile. This trims nucleotides 
        from the sequence being added to the tile until they translate in the 
        same frame.

        Parameters:
        ___________
        hsp : Hsp object if Hsp to add.
        location : either upstream (0) or downstream (1). The location of the 
            Match relative to the existing tile.
        distance : distance between hsp and this tile, used as fill value
        """
        difference = distance
        while difference % 3 != 0: 
            # remove nts until difference is divisible by 3
            difference -= 1

        # add hsp with the cropped nts and translate
        self.add_hsp(hsp, location, difference) 

        print >> logfile, "Hsp number {} not in same frame as tile - correcting".format(hsp.num)

    def extendReadingFrame(self, conservative=False, prokaryotic=False):
        """Expand corrected nucletode sequence to the first start and stop codons
        found in teh same frame

        Parameters:
        ___________
        conservative : bool. Only find first start codon, do not exand to the next 
                       stop farther upstream.
        prokaryotic : bool. Use alternate stop codons found in prokaryotic genomes
        """
        
        if not prokaryotic:
            startCodons = ["ATG"]
        else:
            startCodons = ["ATG", "GTG", "TTG"]

        stopCodons = ["TAG", "TAA", "TGA"]

        #Expand upstream to nearest stop codon (if not conserved) of nearest start codon
        #if conservative.
        codons = startCodons if conservative else stopCodons
        start = self.start
        while start >= 0 and not self.contig.sequence[start:start+3] in codons:
            if start<3:
                break 
            start -= 3

        #Make sure that start is actaully a start or stop codon
        if not self.contig.sequence[start:start+3] in codons:
            start = self.start
            print >> logfile, "Warning, contig {} has no first stop codon within frame".format(self.contig.name)


        #Find closest stop codon
        end = self.end-3
        while end<=len(self.contig.sequence)-3 and not self.contig.sequence[end:end+3] in stopCodons:
            end += 3

        #Make sure that stop is actually a stop codon
        if not self.contig.sequence[end:end+3] in stopCodons:
            end = self.end-3
            print >> logfile, "Warning, contig {} has no stop codon within frame".format(self.contig.name)

        self.start = start
        self.end = end+3

        self.nt_seq.sequence = "{}{}{}".format(self.contig.sequence[start:self.start],
                                               self.nt_seq.sequence,
                                               self.contig.sequence[self.end:end+3])
        self.aa_seq = self.nt_seq.translate_sequence(strand=self.strand)

    def determineGaps(self):
        """Replace Xs with a cartesian product of all nucleotides. The length 
        of Xs signifies the length of the gap. The sequence whos 3-mer count 
        most accuraltey resembles the codon usage table is returned.

        e.g. AAAAXXXXCCGXXCCX

        WARNING: NOT FULLY IMPLEMENTED
        """
        #Store the best sequence
        bestScore = 0
        bestSequence = None

        print "There are {} combinations to try for {}".format(4**self.nt_seq.count("X"), self.contig.name)
        print self.nt_seq

        #Create seval terators that return the position of each X
        xPat = re.compile("X")

        #Try every single combination of nts and see which most closely relates
        #to the codon usage
        for i, replacement in enumerate(product("ACGT", repeat=self.nt_seq.count("X"))):
            xIter = iter(replacement)
            replSeq = xPat.sub(lambda m:xIter.next(), self.nt_seq)

            #Compare this with the codon usage table
            kmers, numSeqs = count_kmers(">t\n{}".format(replSeq).split("\n"), IUPAC_N, 3, normalize=True)
            score = 0
            for kmer, newCount in kmers.iteritems():
                originalUsage = self.codon_usage[kmer]
                score += newCount - originalUsage

            if score < bestScore:
                bestSequence = replSeq


        """Below is another test to add a gap based on codon usage. It first sees if the start
        is a multiple of three and has the length of a multiple of three. If so, try all possible
        codons that fit. If the start is not a multiple of 3, it finds the most probable codons 
        that start with the one or two nucleotides before it. Next is sees if the rest of the gap 
        until the end is a multiple of three. If not, it is clipped until it is. This sequnce is then
        tested with all of the best codons. If the end was not a multiple of three, find codons with 
        one or two nucleotides that follow it. Untested.
        """
        if False:
            for match in re.finditer("X+", self.nt_seq):
                if match.start() % 3 == 0 and len(match.group(0)) % 3 == 0:
                    #Match starts in frame and the length can hold at least one codon
                    for codons in product(self.codon_usage.keys(), repeat=len(match.group(0))/3):
                        replSeq = "{}{}{}".format(self.nt_seq[:match.start()],
                                                  "".join(codons),
                                                  self.nt_seq[match.end():])
                        #Which is the best?
                else:
                    start = match.start()
                    while start % 3 != 0:
                        start -= 1

                    startCodons = [codon for codon in self.codon_usage if codon.startswith(self.nt_seq[start:3].replace("X", ""))]
                    #Which is the best?

                    end = match.end()
                    while end+1 % 3 != 0:
                        end -= 1

                    middle = end-(match.start()-start)+1
                    for codons in product(self.codon_usage.keys(), repeat=middle):
                        replSeq = "{}{}{}".format(self.nt_seq[:match.start()],
                                                  "".join(codons),
                                                  self.nt_seq[match.end():])

                    endCodons = [codon for codon in self.codon_usage if codon.endswith(self.nt_seq[end+1:3].replace("X", ""))]

        self.nt_seq = bestSequence

    def _getCorrectedAA(self, query_seq):
        """Get the aa sequence used as the query and correct sequences that have
        repeptative or low-complexity regions filtered from BLASTX by Removing Xs 
        and gaps.

        This function creates regular expression from the Match query sequence. 
        E.g., if the Match protein sequence were FGTPPXPYII, the regular expression
        created would be FGT....YII. This is used to search all the translations of 
        the contig for a matching sequence.

        Parameters:
        ___________
        query_seq : the query sequence

        Return:
        ________
        aa : corrected aa sequence
        """
        #remove gaps
        aa = re.sub(r'-', r'', query_seq)

        #first 3 aas - replace X with . for re matching  
        start = re.sub(r'[X|x]', r'.', aa[:3]) 

        #change * to _ to match stop codons correctly
        start = re.sub(r'\*', r'.', start)

        #last 3 aas
        end = re.sub(r'[X|x]', r'.', aa[-3:]) 
        end = re.sub(r'\*', r'.', end)

        #Sequences that match one of the sequces in the six frame trnslation
        matches = []
        attempt = 0

        while not matches and attempt < 2:
            #regex for matching
            anychar_pattern = "."*(len(aa)-6)
            pattern = r"{}{}{}".format(start, anychar_pattern, end)

            for frame, protein in self.contig.sixframe():
                #search all protein translations for regex
                match = re.search(pattern, protein.sequence)
                if match:
                    matches.append(match)

            if not matches:
                #if there are no matches, try looking for unknown bases as well
                start= "({}|J)({}|J)({}|J)".format(*start)
                end = "({}|J)({}|J)({}|J)".format(*end)
                attempt += 1

            if attempt == 2:
                raise RuntimeError("Cannot match HSP to Contig {}".format(self.contig.name))

        # and assuming there is only one match, the index is correct
        match = matches[0].group() 
        if len(match) != len(aa): 
            # sanity check
            raise RuntimeError("Cannot match HSP to Contig {}".format(self.contig.name))
        else:
            while 'X' in aa:
                #search for X or runs of X
                xs = re.search(r'X+', aa) 
                
                #indices of xs matched above
                s = xs.start() 
                e = xs.end()

                #portion of contig to use to replace Xs
                subst = match[s:e] 
                
                #only replace first instance in each loop
                aa = re.sub(xs.group(), subst, aa, count=1)
        
        #replace J (from unknown translations) with X as these are genuine unknowns
        aa = re.sub(r'J', r'X', aa)

        return aa

    def printer(self):
        """Writes hsp information to global log File
        """
        print >> logfile, "Sequence: {}".format(self.nt_seq)
        print >> logfile, "Start: {}, End: {}".format(self.start, self.end)
        print >> logfile, "Protein: {}".format(self.aa_seq)

class EmptyTilePath(Tile_Path):
    def __init__(self, contig, protein=False, no_hit=0):
        """Initialise an empty tile path - there were no annotations for this
        contig.

        Parameters:
        ___________
        contig : DNASequence object
        protein : bool. output amino acid sequence
        no_hit : How to hande sequences with no annotation
            0=Ignore
            1=Output longest ORF
            2=Output frame 1
        """
        self.contig = contig
        self.tile = False
        self.start = 0
        self.strand = 1
        self.frame = 1
        self.hitID = None
        self.bitscore = 0.0
        self.start = 0
        self.end = len(self.contig.sequence)
        self.nt_seq = self.contig
        self.aa_seq = self.contig.translate_sequence(strand=self.strand)
        self.saveProteinCoding = True
        self.protein = protein
        self.no_hit = no_hit

    def useLongestORF(self):
        """Find longest ORF in either frame.
        """
        stopCodons = ["TAG", "TAA", "TGA"]

        lenLongestORF = 0
        startLongestORF = 0
        endLongestORF = 0
        frameLongestORF = 1

        openOrf = [False for _ in xrange(3)]

        rc = self.contig.revcomp()

        for strand, seq in zip((+1, -1), (self.contig, rc)):
            for i in xrange(len(seq)):
                frame = i%3
                codon = seq.sequence[i:i+3]
                
                if codon in stopCodons and not openOrf[frame]:
                    #Begin ORF if no other orf is started in current frame 
                    #save position in current frame
                    openOrf[frame] = i

                elif codon in stopCodons and openOrf[frame]:
                    #Stop open ORFS in current frame if an orf is open
                    #If found orf is greater than longest, reset it to the longest
                    orfLength = i-openOrf[frame]
                    if orfLength > lenLongestORF:
                        lenLongestORF = orfLength
                        startLongestORF = openOrf[frame]
                        endLongestORF = i+3
                        frameLongestORF = (frame+1)*strand

                    #Clear open ORF in current frame
                    openOrf[frame] = False
        
        self.start = startLongestORF+3
        self.end = endLongestORF
        self.frame = frameLongestORF
        if frameLongestORF < 0:
            self.nt_seq.sequence = rc.sequence[self.start:self.end]
        else:
            self.nt_seq.sequence = self.contig.sequence[self.start:self.end]
        self.aa_seq = self.contig.translate_sequence(start=self.start, stop=self.stop)

    def write(self, outfile):
        """Save tile to file

        Parameters:
        ___________
        outfile : file-like object to write tile

        Return:
        ________
        tile append to outfile
        """
        if self.no_hit == 1:
            self.useLongestORF()
        elif self.no_hit == 2:
            #This will be done automatically
            pass
        else:
            #Ignore Hit, default
            return

        if self.protein:
            self.outputProtein()
        
        #Save the tiles region
        print >> outfile, self

class HSPTiler(object):
    """Main corrector here
    """
    def __init__(self, 
                 fasta, 
                 annotation,
                 format, 
                 gap_limit=15, 
                 evalue_cutoff=1e-10, 
                 allHits=False, 
                 taxfilter=None,
                 filterType=0,
                 codon_usage=None,
                 threads=1,):
        self.fasta = fasta
        self.annotation = annotation
        self.format = format 
        self.gap_limit = gap_limit
        self.evalue_cutoff = evalue_cutoff 
        self.allHits = allHits
        self.filter = taxfilter
        self.filterType = filterType
        self.codon_usage = codon_usage
        self.threads = threads

        if self.format in ["blastxml", "rapsearch2"]:
            #Initilise BLAST XML parser
            self.parser = BlastXMLParser(self.annotation, allHits=self.allHits, evalue=self.evalue_cutoff, taxFilter=self.filter, taxFilterType=self.filterType)
        #elif format == "pslx": #Not fully implemented yet
            #Initilise BLAT PSLX parser
            #parser = BlatPSLXParser(annotation, allHits=allHits, evalue=evalue_cutoff, taxFilter=filter, taxFilterType=filterType)        
        else:
            raise RuntimeError("Annotation format is not recognized.")

    def run(self, outfile=None):
        """Run HSP-Tiler and yield tiles generated for each contig

        Parameters:
        ___________
        outfile : file-like object to save correct sequences to. Optional.

        Return:
        _______
        tile : a Tile_Path object of the correct sequences
            can be written to a file, or yielded as a generator
        """
        if self.threads > 1:
            fasta1, fasta2 = tee(self.fasta)
            numSequences = sum(1 for line in fasta1 if line.startswith(">"))
            chunks = numSequences/(threads-1)+1
            contigs = read_fasta(fasta2)
            sequence_info = izip_missing(contigs, self.parser.parse(), key=lambda x: x.name, fillvalue=None)

            pool = Pool()
            corrected_sequences = pool.imap(self.correct_frameshifts, sequence_info, chunks)
        else:
            #Read fasta
            contigs = read_fasta(self.fasta)
            # get list of hsps for highest scoring hit for each contig - go to next contig if no hits
            corrected_sequences = imap(self.correct_frameshifts, izip_missing(contigs, self.parser.parse(), key=lambda x: x.name, fillvalue=None))

        if outfile is not None:
            #Write corrected sequences
            for tile in corrected_sequences:
                tile.write(outfile)
        else:
            return corrected_sequences
        

    def correct_frameshifts(self, *args):
        """Correct frameshifts

        Parameters:
        ___________
        contig - DNASequence object to correct
        matches - Matches object

        Returns:
        ________
        tile - corrected sequence
        """
        if len(args) == 1:
            contig, matches = args[0][0], args[0][1]
        elif len(args) == 2:
            contig, matches = args[0], args[1]
        else:
            raise RuntimeError("Wrong number of parameters.")

        print contig.name

        print >> logfile, "\nNext contig:"
        print >> logfile, str(contig)

        #Add each hsp (Match) to list
        hsp_list = list(matches.matches) if matches else None

        if not hsp_list or hsp_list is None:
            print >> logfile, "No hits for {}".format(contig.name)
            contig.description += " [GI=None;S=0.0]"
            return EmptyTilePath(contig)
             
        # Initialise tile using highest scoring hsp (first in list) 
        tile = Tile_Path (hsp_list[0], contig, codon_usage=self.codon_usage) 
        hsp_list[0].used = True # mark that this hsp has been used
        print >> logfile, "Starting tile..."
        hsp_list[0].printer(logfile)
        tile.printer()

        # go to next contig if only one hsp in list
        if len(hsp_list) == 1:
            print >> logfile, "Only one hsp for this hit.  Tile is complete"
            return tile

        # if >1 hsp, attempt to create a tile
        else: 
            # run loop until no hsps were added in the previous iteration
            while tile.added == True: 
                tile.added = False # loop will exit if this is not set to True by calling add.hsp() in this iteration

                # iterate hsps for this contig
                for hsp in hsp_list:

                    # ignore hsps that are already in the tile
                    if hsp.used == False: 
                        print >> logfile, "\nNext hsp..."
                        hsp.printer(logfile)

                        # ignore hsps that are not on the same strand as the tile
                        if hsp.strand != tile.strand: 
                            print >> logfile, "Hsp number {} is not on the same strand as the tile - cannot tile".format(hsp.num)
                            hsp.used = True
                            continue # next hsp

                        # ignore hsps that are completely within the tile
                        if tile.start <= hsp.query_start and hsp.query_end <= tile.end: 
                            print >> logfile, "Hsp {} completely within tile - skipped".format(hsp.num)
                            hsp.used = True # mark as used but do not add to tile
                            tile.printer()
                            continue

                        # determine whether incoming hsp is 5' or 3' of tile
                        # determine whether incoming hsp is 5' or 3' of tile
                        if hsp.query_start <= tile.start: # hsp nt seq 5' of tile
                            location = UPSTREAM
                            distance = abs(hsp.query_end - tile.start)
                        elif hsp.query_end >= tile.end: # hsp nt seq 3' of tile
                            location = DOWNSTREAM
                            distance = abs(hsp.query_start - tile.end)

                        #Attempt to add overlapping sequences that have the same start position as the tile
                        #but is longer than the current start, and adds hsps that have the same end position
                        #but extend longer on the 5' end.
                        """if hsp.query_start < tile.start: # hsp nt seq 5' of tile
                            location = UPSTREAM
                            distance = abs(hsp.query_end - tile.start)
                        elif hsp.query_start == tile.start and hsp.query_end < tile.hsps[0].query_end:
                            location = UPSTREAM
                            distance = abs(hsp.query_end - tile.start)
                        elif hsp.query_start == tile.start and hsp.query_end > tile.hsps[0].query_end:
                            location = DOWNSTREAM
                            distance = abs(hsp.query_start - tile.end)
                        elif hsp.query_end > tile.end: # hsp nt seq 3' of tile
                            location = DOWNSTREAM
                            distance = abs(hsp.query_start - tile.end)
                        elif hsp.query_end == tile.end and hsp.query_start < tile.hsps[-1].query_start:
                            location = UPSTREAM
                            distance = abs(hsp.query_end - tile.start)
                        elif hsp.query_end == tile.end and hsp.query_start > tile.hsps[-1].query_start: # hsp nt seq 3' of tile
                            location = DOWNSTREAM
                            distance = abs(hsp.query_start - tile.end)"""

                        # determine whether hsps overlap and add if they do
                        if ((hsp.query_start <= tile.start <= hsp.query_end) or 
                           (hsp.query_start <= tile.end <= hsp.query_end)):
                            tile.add_overlap(hsp, location, distance) 

                        # otherwise they do not overlap
                        else: 

                            # check hsp locations make sense before continuing
                            if location == UPSTREAM: 
                                assert hsp.query_end < tile.start 
                            elif location == DOWNSTREAM:
                                assert hsp.query_start > tile.end 

                            # check if hsp is close enough to tile to incorporate
                            if distance > self.gap_limit: 
                                print >> logfile, "Too far away to tile"
                                # go to next hsp but do not mark as incorporated as may be able 
                                # to add in a later iteration    
                                continue 

                            # if contig is close enough to attempt tiling                
                            else: 
                                print >> logfile, "Hsp {} is {} nucleotides {} of tile".format(hsp.num, distance-1, location)
                                assert distance -1 >=0, 'Cannot tile hsps' #sanity check

                                #Check whether the gap between the hsps relative to the HIT is 0.  
                                #If this is true but there is a nucleotide gap relative to the contig, 
                                #something must have been inserted into the contig.  Statement is unwieldy 
                                #because test must be made for every combination of hsps upstream and 
                                #downstream of the tile and on positive and negative strands.
                                if ((location == UPSTREAM and 
                                    ((tile.strand == 1 and tile.hsps[0].hit_start - hsp.hit_end -1 == 0) or 
                                     (tile.strand == -1 and hsp.hit_start - tile.hsps[0].hit_end -1 == 0))) or 
                                   (location == DOWNSTREAM and 
                                    ((tile.strand == 1 and hsp.hit_start - tile.hsps[-1].hit_end -1 == 0) or 
                                     (tile.strand == -1 and tile.hsps[0].hit_start - hsp.hit_end -1 == 0)))): 
                                    tile.add_insert(hsp, location, distance)

                                #If the above is not true, check whether the gap between the hsps 
                                #relative to the hit is > 0.  If this is true and there is no nucleotide 
                                #gap relative to the contig, something must have been deleted from the contig.  
                                #If this is true and there is a nucleotide gap, the nucleotides in that gap 
                                #must have been substituted in the contig. Statement is unwieldy as above
                                elif ((location == UPSTREAM and 
                                      ((tile.strand == 1 and tile.hsps[0].hit_start - hsp.hit_end -1 > 0) or 
                                       (tile.strand == -1 and hsp.hit_start - tile.hsps[0].hit_end -1 > 0))) or 
                                     (location == DOWNSTREAM and 
                                      ((tile.strand == 1 and hsp.hit_start - tile.hsps[-1].hit_end -1 > 0) or 
                                       (tile.strand == -1 and tile.hsps[0].hit_end - hsp.hit_start -1 > 0)))): 
                                    tile.add_subst(hsp, location, distance)
                               
                                #as we have already dealt with hsps that overlap the tile, if neither 
                                #of the above are true, there is a problem somewhere
                                else:
                                    print >> logfile, "Contig {}: HSP{} could not be tiled".format(contig.name, hsp.num)
                                    hsp.used = True
                                    continue

        #Fix gaps if codon usage table supplied
        if self.codon_usage and "X" in tile.nt_seq:
            tile.determineGaps()

        # add unchanged contig to output if no changes made. Add tile sequence if changes have been made.
        return tile

def izip_missing(iterA, iterB, **kwds):
    """Iterate through two iterables, while making sure they are in the same
    order. If there are missing values, you can skip the value entirely or 
    return only the iterator with the value and a special fill value.

    Parameters:
    ___________
    iterA : the first iterator
    iterB : the second iterator
    key : function that returns items to compare. Must return strings, ints, or
        an object with the __lt__, __gt__, and __eq__ methods. Optional.
    fillvalue : The value to return if the item is missing. Optional.

    Returns:
    ________
    A : item from first iterator, or fillValue
    B : item from second iterator, or fillValue
    """
    #Get the comparison function
    key = kwds.get("key")
    if key is None:
        key = lambda x: x

    useMissing = False
    fillValue = ""
    if "fillvalue" in kwds:
        useMissing = True
        fillvalue = kwds["fillvalue"]

    #Start both iterators
    A = iterA.next()
    B = iterB.next()

    try:
        while True:
            if key(A) == key(B):
                yield A, B
                A = iterA.next()
                B = iterB.next()
            elif key(A) < key(B):
                if useMissing:
                    yield A, fillvalue
                A = iterA.next()
            elif key(A) > key(B):
                if useMissing:
                    yield fillvalue, B
                B = iterB.next()
            else:
                raise RuntimeError("Invalid compartor")
    except StopIteration:
        pass

#######################################
#Main
#######################################

def parse_args(args):
    """Parsing command line options
    """
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
    parser.add_argument("--format", 
                        default="blastxml", 
                        choices=("blastxml", "blasttab", "pslx", "rs2"),
                        help="Format of the annotation file. Can use output from BLASTX (xml (default), or tabs with specified format), BLAT PSLX, and RapSearch2.")
    # gap limit
    parser.add_argument("-g", "--gap_limit", 
                        type=int, 
                        default=15, 
                        help="Cutoff for distance between hsps. If the gap between hsps in nucleotides is greater than this value, the hsp will not be added to the tile.  Default = 15nt")
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
                        help="Filter out blastx hits from a given species and fmily name using full scientific names only for now")
    parser.add_argument("--filterType",
                        default=0, #0 for ncbi, 1 for uniprot
                        help="Regex to retreive taxon info from hit. Use 0 for NCBI, or 1 for uniprot, or the full regex query")
    #Use threads if there are a lot of sequences
    parser.add_argument("-t", "--threads",
                        default=1,
                        type=int,
                        help="Number of threads to use. Default is 1")
    #Use codon usage to fill gaps?
    codonUsageGroup = parser.add_mutually_exclusive_group()
    codonUsageGroup.add_argument("--codon_usage",
                                help="Path to file, URL, or ID of file contaning codon usage")
    codonUsageGroup.add_argument("--compute_codon",
                                 default=False,
                                 action="store_true",
                                 help="Compute the codon usage table for given fasta sequences")
    #Handle sequences with no hits
    parser.add_argument("--no_hit",
                       choices=range(3),
                       default=0,
                       help="hande sequences with no annotation: 0=Ignore; 1=Output longest ORF; 2=Output frame 1")
    #Define output
    parser.add_argument("-p", "--protein",
                        default=False,
                        action="store_true",
                        help="Output protein sequence. Default is false.")
    parser.add_argument("-o", "--outfile",
                        type=argparse.FileType('wt'),
                        default=sys.stdout,
                        help="File to save corrected sequences")
    parser.add_argument("-l", "--logfile",
                        type=argparse.FileType('wt'),
                        default=sys.stderr,
                        help="File to save log")

    if args[0] == __file__:
        args = args[1:]

    # print help message if no arguments are given
    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    #Parse args
    args = parser.parse_args(args)

    #Create two iterators from same file, one for codon usage other for hsp-tiler
    codon_fasta, args.contigs = tee(args.fasta)

    global logfile
    logfile = args.logfile

    #Parse codon usage
    if args.codon_usage:
        args.codon_usage = read_codon_frequencies(args.codon_usage)
    elif args.compute_codon:
        args.codon_usage = count_kmers(codon_fasta, IUPAC_N, 3, normalize=True)
    else:
        args.codon_usage = None

    return args

def main(args):
    #Parse arguments
    args = parse_args(args)

    #Save file info to log
    print >> logfile, "# hsp_tiler.py output file"
    print >> logfile, "# Date: {}".format(datetime.now())
    print >> logfile, "# Fasta: {}".format(args.fasta.name)
    print >> logfile, "# Annotation: {}".format(args.annotation.name)
    print >> logfile, "# Outfile: {}".format(args.outfile.name)

    if args.codon_usage:
        print >> logfile, "#Codon Usage: {}".format(args.codon_usage)

    print >> logfile, "# Parameters: gap_limit={}, evalue_cutoff={}, useHit1={}, filter=\"{}\"".format(args.gap_limit,
                                                                                                       args.evalue_cutoff,
                                                                                                       args.useHit1,
                                                                                                       args.filter)

    


    hsp_tiler = HSPTiler(args.contigs, 
                         args.annotation,
                         args.format, 
                         gap_limit=args.gap_limit,
                         evalue_cutoff=args.evalue_cutoff, 
                         allHits=not args.useHit1, 
                         taxfilter=args.filter,
                         filterType=args.filterType,
                         codon_usage=args.codon_usage)

    hsp_tiler.run(args.outfile)

    #Cleanup
    logfile.close()
    args.outfile.close()

if __name__ == "__main__":
    main(sys.argv)

# Author: Eli Draizen
# Date 5-1-14
# File: BlatPSLXParser.py

"""Parser for BLAT PSLX file. BLAT is faster than BLASTX for many cases so it
may be an alternative to those who need to correct sequences in a time limit.

To get BLAT to work like BLASTX, you must manually convert the contigs into the
six frames and rerun blast. There is a program in HSP-Tiler utiliities to do
this called blatx.py 

Command:
blatc.py contigs.fa proteins.fa
blat -prot -out=pslx nr.fa proteins.fa
"""

#Standard libraries
import csv
from collections import defaultdict
from itertools import groupby, ifilter, imap

#Custom libraires
from Parser import Parser, Match, Matches

class PSL(Match):
    def __init__(self, *args):
        """Create a PSL object with 23 values in the follwoing order:
        matches - Number of matching bases that aren't repeats.
        misMatches - Number of bases that don't match.
        repMatches - Number of matching bases that are part of repeats.
        nCount - Number of 'N' bases.
        qNumInsert - Number of inserts in query.
        qBaseInsert - Number of bases inserted into query.
        tNumInsert - Number of inserts in target.
        tBaseInsert - Number of bases inserted into target.
        strand - defined as + (forward) or - (reverse) for query strand. In mouse, a second '+' or '-' indecates genomic strand.
        qName - Query sequence name.
        qSize - Query sequence size.
        qStart - Alignment start position in query.
        qEnd - Alignment end position in query.
        tName - Target sequence name.
        tSize - Target sequence size.
        tStart - Alignment start position in query.
        tEnd - Alignment end position in query.
        blockCount - Number of blocks in the alignment.
        blockSizes - Comma-separated list of sizes of each block.
        qStarts - Comma-separated list of start position of each block in query.
        tStarts - Comma-separated list of start position of each block in target.
        qSeq - Query sequence
        tSeq - Target sequence 

        Rember that to Datbase and query and switch for BLAT
        """
        if not len(args) == 23:
            raise RuntimeError("Must include all specidied parameters")

        #Parse PSL
        self.matches = int(args[0])
        self.misMatches = int(args[1])
        self.repMatches = int(args[2])
        self.nCount = int(args[3])
        self.qNumInsert = int(args[4])
        self.qBaseInsert = int(args[5])
        self.tNumInsert = int(args[6])
        self.tBaseInsert = int(args[7])
        self.qName = args[9]
        self.qSize = args[10]
        self.qStart = int(args[11])
        self.qEnd = int(args[12])
        self.tName = args[13]
        self.tSize = int(args[14])
        self.tStart = int(args[15])
        self.tEnd = int(args[16])
        self.blockCount = int(args[17])
        self.blockSizes = map(int, args[18].split(",")[:-1])
        self.qStarts = map(int, args[19].split(",")[:-1])
        self.tStarts = map(int, args[20].split(",")[:-1])
        self.qSeq = args[21]
        self.tSeq = args[22]

        #Rename variables for HSP-Tiler
        self.strand = 1 if args[8][1] == '+' else -1
        self.frame = self.getFrame()
        self.query_start = self.tStart-1 #position in nts relative to query => the db tStart
        self.query_end = self.tEnd #tEnd
        self.hit_start = self.qStart-1 #position in aas relative to hit -> the query qStart
        self.hit_end = self.qEnd #qEnd
        self.bitscore = self.getScore()
        self.used = False # True marks hsps that have already been incorporated into the tile
        self.num = 0
        self.query_seq = self.tSeq

        try:
            #Assume NCBI headers
            self.hitID = self.qName.split("|")[1] #field is 'gi|XXXXX|ref|abcd', save XXXX
        except:
            self.hitID = self.qName #just use the whole hit id

    def isProtein(self):
        """This function automatically determines whether or not the PSL output 
        file contains alignment information for a protein query. 
        is psl a protein psl (are it's blockSizes and scores in protein space)
        Converted from Kent's C code into python

        Parameters:
        ___________
        psl : PSL object

        Return:
        _______
        True if psl is protein
        """

        if self.strand > 0:
            return self.tEnd == self.tStarts[-1] + 3*self.blockSizes[-1]
        else:
            return self.tStart == (self.tSize-(self.tStarts[-1] + \
                3*self.blockSizes[-1]))

    def getScore(self):
        """Return score for psl
        Converted from Kent's C code into python
        """
        sizeMul = 3 if self.isProtein() else 1

        return sizeMul * (self.matches + (self.repMatches>>1)) - \
         sizeMul * self.misMatches - self.qNumInsert - self.tNumInsert

    def getFrame(self):
        """
        """
        return (self.tStart%3+1)*self.strand

class BlatPSLXParser(Parser):
    """Parses BLAT PSLX files, and yields Match objects fora given contig
    """
    def __init__(self, infile, contigs, *args, **kwds):
        """Initialize PSLX parser. Save the contigs to get the order of PSLs

        Parameters:
        ___________
        infile : file-like object containg BLAST XML output
        contigs : list of contigs
        """
        Parser.__init__(self, infile, **kwds)
        self.order = map(lambda f: f.name, contigs)

    def _parsePSLX(self):
        """Read in PSLX file, yeilding PSL objects
        """

        for line, fields in enumerate(csv.reader(self.infile, delimiter="\t")):
            if line < 5:
                continue
            
            yield PSL(*fields)

    def parse(self):
        """Real work done here, converts PSLs into Matches and return them.
        Also switches query and target to be parallel with BLASTX and can used
        with HSP-Tiler.

        Yields:
        matches - Matches object with contig name and Match objects (HSP equivalent)
        """
        print "Reading PSLX..."
        for psl in groupby(self._parsePSLX(), key=lambda psl: psl.qName):
            if self.taxFilter and self._filterTaxa(psl.qName):
                continue

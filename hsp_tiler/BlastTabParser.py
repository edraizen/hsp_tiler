# Author: Eli Draizen
# Date 5-1-14
# File: BlastXMLParser.py

#Standard libraries
import sys
from math import exp
import csv
from itertools import imap, groupby

#Custom libraries
from Parser import Parser, Match, Matches

class HSP(Match):
    """Holds information about a given High-scoring Sequence Pair (HSP) 
    returned from BLASTX.
    """
    def __init__(self, hsp, qframe=True, query_seq=True):
        """Initialise an Hsp object, with the XML 
        representation of the hsp, and a contig.  

        Parameter:
        __________
        hsp : list of atttributes in the following order: qseqid, sseqid, pident, 
            length, mismatch, gapopen, qstart, qend, sstart, send, evalue, 
            bitscore, (qframe, score, qseq).
        """
        self.frame = int(hsp[12]) if qframe else None
        self.query_start = int(hsp[6])-1 #position in nts relative to query
        self.query_end = int(hsp[7])
        self.query_length = self.query_end-self.query_start+1
        self.strand = 1 if self.frame > 0 else -1


        self.hit_start = int(hsp[8])-1 #position in aas relative to hit
        self.hit_end = int(hsp[9])
        self.hit_length = self.hit_end-self.hit_start+1

        self.evalue = float(hsp[10])

        self.bitscore = float(hsp[11])

        self.hitID = self._getGI(hsp[1])
        
        self.num = 0

        if query_seq:
            self.query_seq = hsp[12+int(qframe)+int(score)]
        else:
            self.query_seq = ""

        self.used = False # True marks hsps that have already been incorporated into the tile

    def _getGI(self, description):
        try:
            #Assume NCBI headers
            hitID = description.split("|")[1] #field is 'gi|XXXXX|ref|abcd', save XXXX
        except:
            hitID = description #just use the whole hit id

        return hitID

class BlastTabParser(Parser):
    """Parses BLAST XML files and easily separates Iterations, Hits, and HSPs without 
    loading all of the data into memory. Also can filter results by evalue or
    Taxonmic information.

    Reinventing the wheel. Biopython may be more robust, but this does the job.
    """
    def __init__(self, *args, **kwds): #blast, allHits=False, evalue=1e-10, taxFilter=None, taxFilterType=0):
        """Intialise a BLAST XML parser.
        """
        #Boolean to allow hits to be returned
        self.runHSP = True
        #Save the current contig, or queryID 
        self.queryID = None
        #The number of Hits that have been processed
        self.numHits = 0

        self.qframe = kwds.get("qframe", True)
        self.score = kwds.get("score", False)
        self.query_seq = kwds.get("query_seq", False)

        if "qframe" in kwds:
            del kwds["qframe"]

        if "score" in kwds:
            del kwds["score"]

        if "query_seq" in kwds:
            del kwds["query_seq"]

        #Initialize Parser super class
        Parser.__init__(self, *args, **kwds)

    def parse(self):
        """Process blastx tab delimeted output. Only saves best hit for each query
        and only save queries who are in the the list of ids

        Input:
        -resultsFile, file-like object of blastx output
        -ids, list of ids of query ID's

        Output:
        -scores, list of bitsocres for each hit
        """
        def skip_comments(iterable):
            for i, line in enumerate(iterable):
                if not line.startswith('#'):
                    yield line

        for contig, hits in groupby(csv.reader(skip_comments(self.infile), delimiter="\t"), key=lambda x:x[0]):
            results = [HSP(x, qframe=self.qframe, query_seq=self.query_seq) for x in hits]
            #map(lambda x: HSP(x, qframe=self.qframe, query_seq=self.query_seq), hits)
            yield Matches(contig, results)


# Author: Eli Draizen
# Date 5-3-14
# File: Parser.py

"""Module to handle Matches from an annotation program such as BLASTX or BLAT

Can be used to build parsers for new programs
"""

#Standard libraries
import re

#Parse out taxa info
regexTaxFilters = {0:re.compile("^.+\[(.+)\]"),  #nr
                   1:re.compile("OS=(\w+\s\s+)") #uniprot
                  }

class Parser(object):
    """Base class for building new parsers. 
    """
    def __init__(self, infile, allHits=True, evalue=1e-10, taxFilter=None, taxFilterType=0):
        """Intialise a BLAST XML parser. 

        Parameters:
        ___________
        infile : file-like object containg BLAST XML output
        allHits : Use all of the hits instead of first, top hit
        evalue : evalue cutoff; defialt 1e-10
        taxFilter : Ignore hits from given species or genus
        taxFilterType : Regex to parse taxon information, or 0 for nr, or 1 for uniprot.
        """
        self.infile = infile
        #Use all of the HSP instead of first, top hit
        self.allHits = allHits
        #Save species name or genus and ignore hits from them
        self.taxFilter = taxFilter
        #E-value cutoff
        self.evalue = evalue

        self.taxFilterType = taxFilterType
        self.taxRegex = None

    def _filterTaxa(self, taxaString):
        """Return true if string does not contain taxon filter, else False
        """
        if self.taxRegex is None:
            if not self.taxFilterType in regexTaxFilters:
                self.taxFilterType = re.compile(taxFilterType)
            else:
                self.taxFilterType = regexTaxFilters[taxFilterType]

        try:
            #Parse out genus and species if exists
            taxInfo = self.taxFilterType.findall(taxaString)[0] 
            if len(taxInfo.split()) == 1: 
                #Error finding taxonomic info
                taxInfo = taxaString
        except:
            taxInfo = taxaString
        finally:
            #Only return HSP if the filter is not found within the taxInfo string
            return self.taxFilter not in taxInfo

    def parse(self):
        """Start parsing file format. Must be subclassed to add functionality.

        Yields:
        matches - Matches object with contig name and Match objects (HSP equivalent)
        """
        raise NotImplementedError


class Match(object):
    """Holds information about a given match. The  information held is the 
    strand of the query that the HSP is on, the frame it is in relative to the 
    query sequence, the start and end positions of the HSP relative to the 
    query sequence and the start and end positions relative to the protein 
    sequence of the hit, the e-value and bitscore of the hit, the amino acid 
    sequence of the HSP and the number of the HSP. This object also keep track 
    to see if HSP has already be tiled.

    Must be subclassed to add required information specified by each format.
    """
    #Required information for the Match, set in init!
    frame = 1
    strand = 1
    query_start = 0
    query_end = 0
    hit_start = 0
    hit_end = 0
    evalue = 0
    bitscore = 0
    hitID = ""
    used = False 
    num = 0
    query_seq = ""

    def __init__(self, hsp):
        """Initialise a Match object. Must override this method to set the
        required Match information above

        Parameters:
        ________
        hsp : Element object containing HSP information
        """
        raise NotImplementedError
        
    def printer (self, logfile):
        """Writes hsp information to global log File

        Parameters:
        ___________
        logfile : file-like object to write log
        """
        print >> logfile, "Hsp {}".format(self.num)
        print >> logfile, "Strand: {}".format(self.strand)
        print >> logfile, "Frame: {}".format(self.frame)
        print >> logfile, "Start: {}, End: {}".format(self.query_start, self.query_end)

class Matches(object):
    """Store matches with the contig name they come from"""
    def __init__(self, name, matches):
        """
        Paramters:
        __________
        name : name of contig
        matches : iterable of Match objects
        """
        self.name = name
        self.matches = matches

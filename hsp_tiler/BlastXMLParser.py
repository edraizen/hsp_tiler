# Author: Eli Draizen
# Date 5-1-14
# File: BlastXMLParser.py

#Standard libraries
import sys
import xml.etree.cElementTree as ET
import re
from math import exp

#Custom libraries
from Parser import Parser, Match, Matches

class HSP(Match):
    """Holds information about a given High-scoring Sequence Pair (HSP) 
    returned from BLASTX.
    """
    def __init__(self, hsp, rapsearch=False):
        """Initialise an Hsp object, with the XML 
        representation of the hsp, and a contig.  

        Parameter:
        __________
        hsp : Element object containing HSP information
        """
        frame = int(hsp.find('Hsp_query-frame').text)
        if rapsearch:
        	#1 = frames 0,1,2; -1 = frames 3,4,5
        	self.frame = (frame%3)+1
        	self.strand = 1 if frame < 3 else -1

        	if frame >= 3:
        		#RAPSearch2 switches start and from, switch back
        		self.query_start = int(hsp.find('Hsp_query-to').text)-1 #position in nts relative to query
        		self.query_end = int(hsp.find('Hsp_query-from').text)
        	else:
        		self.query_start = int(hsp.find('Hsp_query-from').text)-1 #position in nts relative to query
        		self.query_end = int(hsp.find('Hsp_query-to').text)
        else:
	        # 1 = frames 1,2,3; -1 = frames -1,-2,-3 
	        self.frame = abs(frame)
	        self.strand = 1 if 0<frame<=3 else -1
        	self.query_start = int(hsp.find('Hsp_query-from').text)-1 #position in nts relative to query
        	self.query_end = int(hsp.find('Hsp_query-to').text)
        self.query_length = self.query_end-self.query_start+1

        self.hit_start = int(hsp.find('Hsp_hit-from').text)-1 #position in aas relative to hit
        self.hit_end = int(hsp.find('Hsp_hit-to').text)
        self.hit_length = self.hit_end-self.hit_start+1

        try:
        	self.evalue = float(hsp.find('Hsp_evalue').text)
        except:
        	try:
        		#Already converted evalue in parser
        		self.evalue = float(hsp.find('Hsp_log-evalue').text)
        	except:
        		raise RuntimeError("XML file must be from BLASTX or RAPSearch2")

        self.score = float(hsp.find('Hsp_score').text)
        self.bitscore = float(hsp.find('Hsp_bit-score').text)
        self.hitID = hsp.find("hitID").text
        self.used = False # True marks hsps that have already been incorporated into the tile
        self.num = int(hsp.find('Hsp_num').text)
        self.query_seq = hsp.find('Hsp_qseq').text

class BlastXMLParser(Parser):
	"""Parses BLAST XML files and easily separates Iterations, Hits, and HSPs without 
	loading all of the data into memory. Also can filter results by evalue or
	Taxonmic information.

	Reinventing the wheel. Biopython may be more robust, but this does the job.
	"""
	def __init__(self, *args, **kwds): #blast, allHits=False, evalue=1e-10, taxFilter=None, taxFilterType=0):
		"""Intialise a BLAST XML parser.
		"""
		#Initialize Parser super class
		Parser.__init__(self, *args, **kwds)

		#Build iter to loop over XML
		self.context = iter(ET.iterparse(self.infile, events=("start", "end")))
		#Boolean to allow hits to be returned
		self.runHSP = True
		#Save the current contig, or queryID 
		self.queryID = None
		#The number of Hits that have been processed
		self.numHits = 0
		#File came form RAPSearch2?
		self.rapsearch = kwds.get("rapsearch", False)

		#Start initial parsing
		event, root = self.context.next()

		if root.tag not in ["BlastOutput", "Output"]:
			raise RuntimeError("This is not a valid BLAST XML file or RAPSearch2 XML file")
		elif root.tag == "Output":
			self.rapsearch = True

		#Start looping over data until we get to first iteration
		for event, elem in self.context:
			if event == "start" and elem.tag == "Iteration":
				break 

	def parse(self):
		"""Real work done here. Combine the query with its HSPs
		"""
		for query in self.parseQuery():
			yield Matches(query, self.parseHsp())


	def parseQuery(self):
		"""Parse each query (i.e. Iteration and Hits), but only use first hit unless allHits is specified
		Sets up parser for user to select HSP with the parseHSP method.

		Return: Contig name or QueryID of the current query. 
		"""
		for event, elem in self.context:
			if event == "end" and elem.tag == "Iteration_query-def":
				#Save current contig or queryID and 
				self.queryID = elem.text
				yield self.queryID
			if self.runHSP and event == "start" and elem.tag == "Hit":
				self.runHSP = self.allHits #Continue looking for hits if allHits is True
			if not self.runHSP and event == "end" and elem.tag == "Iteration_hits":
				self.runHSP = True

	def parseHsp(self):
		"""Process each HSP for a given hit. Must be called during parseQuery.

		Return: XML object of HSP. If allHits is specified, all HSPS for every
		hit are returned for each contig, otehrwise only the HSP within the first
		hit are used.
		"""

		#Decide whether the current HSP is returned
		returnHsp = False

		#Save ID for hit
		hitID = None

		for event, elem in self.context:
			if event == "end" and elem.tag == "Hit_id":
				hitID = self._getGI(elem.text)
			elif event == "end" and elem.tag == "Hit_def":
				#Get species and genus for use with filter

				#One more hit has been seen
				self.numHits += 1

				if hitID is None:
					#RAPSearch doesn't use Hit_id
					hitID = self._getGI(elem.text)

				#Returns all HSP if there is no fitler, else only return HSPS that do not contain filter
				returnHsp = self._filterTaxa(elem.text) if self.taxFilter else True

			elif event == "end" and elem.tag == "Hsp_evalue" and float(elem.text) > self.evalue:
				#Don't return HSP is evalue is above cutoff
				returnHsp = False
			elif event == "end" and elem.tag == "Hsp_log-evalue":
				#RAPSearch2 uses log(evalue)
				evalue = exp(float(elem.text))
				elem.text = evalue
				if evalue > self.evalue:
					returnHsp = False

			elif returnHsp and event == "end" and elem.tag == "Hsp":
				#Add the hitID to current HSP to process later
				hitElement = ET.SubElement(elem, "hitID")
				hitElement.text = hitID
				hsp = HSP(elem, rapsearch=self.rapsearch)
				yield hsp
			elif not self.allHits and event == "end" and elem.tag == "Hit_hsps":
				#Stop after 1st HSP hits are finished
				break
			elif self.allHits and event == "end" and elem.tag == "Iteration_hits":
				#Stop after all of HSPs from every hit are finished
				break

	def _getGI(self, description):
		try:
			#Assume NCBI headers
			hitID = description.split("|")[1] #field is 'gi|XXXXX|ref|abcd', save XXXX
		except:
			hitID = description #just use the whole hit id

		return hitID

if __name__ == "__main__":
	try:
		blast = open(sys.argv[1])
	except:
		raise RuntimeError("Cannot open BLAST XML file")

	parser = BlastXMLParser(blast)

	for queryID in parser.parseQuery():
		print queryID
		for hsp in parser.parseHsp():
			print "\t", hsp
			print "\t\t", abs(int(hsp.find('Hsp_query-frame').text))
        	print "\t\t", int(hsp.find('Hsp_query-from').text) #position in nts relative to query
        	print "\t\t", int(hsp.find('Hsp_query-to').text)


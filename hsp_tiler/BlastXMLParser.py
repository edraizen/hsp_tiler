import sys
import xml.etree.cElementTree as ET
import re

regexTaxFilters = {0:re.compile("^.+\[(.+)\]"),  #nr
                   2:re.compile("OS=(\w+\s\s+)") #uniprot
                  }

class BlastXMLParser(object):
	"""Parses BLAST XML files and easily separates Iterations, Hits, and HSPs without 
	loading all of the data into memory. Also can filter results by evalue or
	Taxonmic information.

	Reinventing the wheel. Biopython may be more robust, but this does the job.
	"""
	def __init__(self, blast, allHits=False, evalue=1e-10, taxFilter=None, taxFilterType=0):
		"""Intialise a BLAST XML parser. 

		Parameters:
		blast - a file-like object containg BLAST XML output
		allHSP - Use all of the HSP instead of first, top hit
		evalue - evalue cutoff; defialt 1e-10
		taxFilter - Ignore hits from given species or genus
		"""
		#Use all of the HSP instead of first, top hit
		self.allHits = allHits
		#Save species name or genus and ignore hits from them
		self.taxFilter = taxFilter
		#E-value cutoff
		self.evalue = evalue

		#Build iter to loop over XML
		self.context = iter(ET.iterparse(blast, events=("start", "end")))
		#Boolean to allow hits to be returned
		self.runHSP = True
		#Save the current contig, or queryID 
		self.queryID = None
		#The number of Hits that have been processed
		self.numHits = 0

		#Start initial parsing
		event, root = self.context.next()
		if not root.tag == "BlastOutput":
			raise RuntimeError("This is not a valid BLAST XML file")

		#Start looping over data until we get to first iteration
		for event, elem in self.context:
			if event == "start" and elem.tag == "Iteration":
				break

		if taxFilter:
			if not taxFilterType in regexTaxFilters:
				self.taxFilterType = re.compile(taxFilterType)
			else:
				self.taxFilterType = regexTaxFilters[taxFilterType]


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
				hitID = elem.text.split("|")[1] #field is 'gi|XXXXX|ref|abcd', save XXXX
			elif event == "end" and elem.tag == "Hit_def":
				#Get species and genus for use with filter

				#One more hit has been seen
				self.numHits += 1

				#Returns all HSP if there is no fitler, else only return HSPS that do not contain filter
				returnHsp = not self.taxFilter
				if self.taxFilter:
					try:
						#Parse out genus and species if exists
						taxInfo = self.taxFilterType.findall(elem.text)[0]
						if len(taxInfo.split()) == 1: 
							#Error finding taxonomic info
							taxInfo = elem.text
					except:
						taxInfo = elem.text
					finally:
						#Only return HSP is the filter is not found within the taxInfo string
						returnHsp = self.taxFilter not in taxInfo

			elif event == "end" and elem.tag == "Hsp_evalue" and float(elem.text) > self.evalue:
				#Don't return HSP is evalue is above cutoff
				returnHsp = False

			elif returnHsp and event == "end" and elem.tag == "Hsp":
				hitElement = ET.SubElement(elem, "hitID")
				hitElement.text = hitID
				yield elem
			elif not self.allHits and event == "end" and elem.tag == "Hit_hsps":
				#Stop after 1st HSP hits are finished
				break
			elif self.allHits and event == "end" and elem.tag == "Iteration_hits":
				#Stop after all of HSPs from every hit are finished
				break

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



import sys
import xml.etree.cElementTree as ET
import re

class BlastXMLParser(object):
	def __init__(self, blast, allHits=False, taxFilter=None):
		"""Intialise a BLAST XML parser. Easily separate Iterations, Hits, and
		HSPs without loading all of the data into memory.

		Parameters:
		blast - a file-like object containg BLAST XML output
		allHSP - Use all of the HSP instead of first, top hit
		filter - Ignore hits from given species or genus
		"""
		#Use all of the HSP instead of first, top hit
		self.allHits = allHits
		#Save species name or genus and ignore hits from them
		self.taxFilter = taxFilter
		#Build iter to loop over XML
		self.context = iter(ET.iterparse(blast, events=("start", "end")))
		#Boolean to allow hits to be returned
		self.runHSP = True

		#Start initial parsing
		event, root = self.context.next()
		if not root.tag == "BlastOutput":
			raise RuntimeError("This is not a valid BLAST XML file")

		#Start looping over data until we get to first iteration
		for event, elem in self.context:
			if event == "start" and elem.tag == "Iteration":
				break

	def parseQuery(self):
		"""Parse each query (i.e. Iteration and Hits), but only use first hit
		"""
		for event, elem in self.context:
			if event == "end" and elem.tag == "Iteration_query-def":
				queryID = elem.text
				yield queryID
			if self.runHSP and event == "start" and elem.tag == "Hit":
				self.runHSP = self.allHits #Continue looking for hits if allHits is True
			if not self.runHSP and event == "end" and elem.tag == "Iteration_hits":
				self.runHSP = True

	def parseHsp(self):
		"""Process each HSP for a given hit. Must be called during parseQuery.
		Return: XML object of HSP
		"""
		returnHsp = False
		for event, elem in self.context:
			if event == "end" and elem.tag == "Hit_def":
				returnHsp = not self.taxFilter
				if self.taxFilter:
					try:
						taxInfo = re.findall("^.+\[(.+)\]", elem.text)[0]
						if self.taxFilter not in taxInfo:
							returnHsp = True
					except:
						pass

			if returnHsp and event == "end" and elem.tag == "Hsp":
				yield elem
			if event == "end" and elem.tag == "Hit_hsps":
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


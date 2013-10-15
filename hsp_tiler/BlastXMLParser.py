
import sys
import xml.etree.cElementTree as ET

class BlastXMLParser(object):
	def __init__(self, blast):
		"""Intialise a BLAST XML parser. Easily separate Iterations, Hits, and
		HSPs without loading all of the data into memory.

		Parameters:
		blast - a file-like object containg BLAST XML output
		"""
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
				continue
			if self.runHSP and event == "start" and elem.tag == "Hit":
				self.runHSP = False
				yield queryID
			if not self.runHSP and event == "end" and elem.tag == "Iteration_hits":
				self.runHSP = True

	def parseHsp(self):
		"""Process each HSP for a given hit. Must be called during parseQuery.
		"""
		for event, elem in self.context:
			if event == "end" and elem.tag == "Hsp":
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


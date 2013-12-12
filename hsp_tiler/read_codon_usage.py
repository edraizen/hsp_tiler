#!/usr/bin/env python2.7
#Author: Eli Draizen
#BME 205 Assignment 10
#File: read_codon_usage

"""Module to read in a file containing codon usage in the format specified by
http://www.kazusa.or.jp/'s Codon Usage Table with Amino Acids format.
"""

import re
import string
import urllib2
from collections import defaultdict

def read_codon_usage(usageFile):
    """read in a file containing codon usage in the format specified by
http://www.kazusa.or.jp/'s Codon Usage Table with Amino Acids format. Returns a
dictionary with codons as a keys and a dictionary of the amino acid it codes for,
the fraction, freq, and number that it occurs.

    usageFile - can be a file-like object of or path/url to a codon usage table
    
    fields: [triplet] [amino acid] [fraction] [frequency: per thousand] ([number])
    Example input:
    UUU F 0.58 22.2 ( 35846)  UCU S 0.14  8.7 ( 14013)  UAU Y 0.57 16.5 ( 26648)
    """
    if not hasattr(usageFile, "read"):
        #Not a file-like object
        try:
            usageFile = open(usageFile)
        except:
            print "here"
            try:
                usageFile = read_codon_web(usageFile)
            except:
                RuntimeError("Unable to load Codon Usage File {}".format(usageFile))
    
    pat = re.compile("""([A-Z]{3})            #Match codon
                         \s+        
                         ([A-Z*]{1})           #Match amino acid
                         \s+
                         ([0-9]+(?:\.[0-9]*)?) #Match fraction
                         \s+ 
                         ([0-9]+(?:\.[0-9]*)?) #Match frequency
                         \s*
                         \(
                            (.+?)           #Match number
                        \)""", re.VERBOSE)
    u_to_t = string.maketrans("U", "T")
    usage = {}
    for line in usageFile:
        line = line.rstrip()
        print line
        for codon, aa, frac, freq, number in pat.findall(line):
            usage[codon.translate(u_to_t)] = {"aa": aa, 
                                                "frac": float(frac), 
                                                "freq": float(freq), 
                                                "number": float(number)}
    return usage

def read_codon_web(url):
    if str(url).isdigit():
        url = "http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species={}\
&aa=1&style=N".format(url)
    elif "kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=" in url:
        url = "{}&aa=1&style=N".format(url)
    else:
        raise RuntimeError("Error! URL not recognized. Must be from kazusa.or.jp/\
codon/cgi-bin/showcodon.cgi")
    codon_usage_scrape = urllib2.urlopen(url)
    codon_usage = ""
    reading_codon_usage = False

    for line in codon_usage_scrape:
        if line.startswith("<PRE>"):
            reading_codon_usage = True
        elif line.startswith("</PRE>"):
            break
        elif reading_codon_usage:
            codon_usage += line
    
    return codon_usage.splitlines()

def read_codon_frequencies(usageFile):
    codon_bias = read_codon_usage(usageFile)
    print codon_bias.keys()
    totalCodonCount = sum([info["number"] for codon, info in codon_bias.iteritems()])
    t = defaultdict(lambda: 1.0/totalCodonCount)
    for codon, info in codon_bias.iteritems():
        t[codon] = info["number"]/float(totalCodonCount)
    print t.keys()
    return t
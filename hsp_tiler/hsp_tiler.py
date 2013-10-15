#!/usr/local/bin/python3
#Author: Kathryn Crouch

#######################################
#Imports
import xml.etree.ElementTree as ET
import re
from collections import deque
import argparse
import sys
#######################################
#Classes

#Contig class
class Contig(object):
    #Initialise - contig is a single fasta item 
    def __init__(self, contig):
        self.contig = contig
        self.name = self.name()
        self.sequence = self.sequence()
        self.length = self.length()
        self.sixframe = self.sixframe()

    # returns contig header text
    def name (self): 
        name = self.contig.split('\n',1)[0].rstrip('\n').lstrip()#split at first newline
        return name

    # returns contig sequence
    def sequence (self): 
        sequence = self.contig.split('\n',1)[1].rstrip('\n')
        return sequence

    # returns sequence length
    def length (self): 
        length = len(self.sequence)
        return length

    # returns a dict containing contig translation in 6 frames
    def sixframe (self): 
        rev_comp = revcomp(self.sequence)
        sixframe = {}
        for i in range (0, 3):
            seq = self.sequence[i:] # get frames positive strand
            protein = codons(seq)
            sixframe [i+1] = protein
            rev = rev_comp[i:] # reverse complement and do again to get the other 3 frames
            protein = codons(rev)
            sixframe [0-i-1] = protein
        return sixframe

    # writes contig information
    def printer (self): 
        logfile.write(self.name + '\n')
        logfile.write (self.sequence + '\n')
        logfile.write ('Length: ' + str(self.length) + '\n')

#######################################        

#XML Handling
class Blast_Handler(object):
    #Initialise - tree is an ElementTree object 
    def __init__(self, tree):
        self.tree = tree         
        self.parent_map = self.parent_child()
        self.queries = list(self.query_list()) #turn iterator into a list

    # creates dict in the form {child:parent} for all nodes in tree
    def parent_child (self): 
        parent_map = dict((c, p) for p in self.tree.iter() for c in p)
        return parent_map

    # returns list of all queries in xml file by contig name
    def query_list (self): 
        queries = self.tree.iter('Iteration_query-def')
        return queries

    # returns highest scoring (first) hit for query specified by contig name
    def get_hit (self, name): 
        for query in self.queries:
            if query.text == name:
                query_parent = self.parent_map[query] # Find parent because Iteration_query-def has no children
                hits = query_parent.findall('.//Hit') # Use xpath - findall only finds direct children
                if hits:
                    hit = hits[0]
                else:
                    hit = None
        return hit

    # returns a list of hsps for the first hit of the query specified by the contig name
    def get_hsps (self, name): 
        hit = self.get_hit(name)
        if hit != None: 
            self.hsps = hit.findall('.//Hsp')
        else:
            self.hsps = None
        return self.hsps

#######################################

#HSPs
class Hsp(object):
    #Initialise from xml - could write another constructor to initialise from another format
    def __init__(self, hsp, contig):
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
        self.num = int(hsp.find('Hsp_num').text)
        self.aa_seq = self.remove_x(hsp,contig) 
        self.used = False # True marks hsps that have already been incorporated into the tile
        

    # Removes Xs and gaps from hsp aa sequences by substituting from contig translation
    def remove_x (self, hsp, contig): 
        no_x = re.sub(r'-', r'', hsp.find('Hsp_qseq').text) # remove gaps 
        start = re.sub(r'[X|x]', r'.', no_x[:3]) # first 3 aas - replace X with . for re matching
        start = re.sub(r'\*', r'.', start) # change * to _ to match stop codons correctly
        end = re.sub(r'[X|x]', r'.', no_x[-3:]) # last 3 aas
        end = re.sub(r'\*', r'.', end) 
        print('Trying to match HSP ' + str(self.num) + ' to Contig ' + contig.name)
        print(no_x)
        for sequence in contig.sixframe.keys():
            print(contig.sixframe[sequence])
        matches = self.match(no_x, contig, start, end)
        if len(matches) != 1: # if there are no matches, try looking for unknown bases as well
            start='(' + start[0] + '|J)(' + start[1] + '|J)(' + start[2] + '|J)'
            end = '(' + end[0] + '|J)(' + end[1] + '|J)(' + end[2] + '|J)'
            matches = self.match(no_x, contig, start, end)
        match = matches[0].group() # and assuming there is only one match, the index is correct
        if len(match) != len(no_x): # sanity check
            raise SystemExit ('Cannot match HSP ' + str(hsp.num) + ' to Contig ' + contig.name)
        else:
            while 'X' in no_x:
                xs = re.search(r'X+', no_x) # search for X or runs of X
                s = xs.start() # indices of xs matched above
                e = xs.end()
                subst = match[s:e] # portion of contig to use to replace Xs
                no_x = re.sub(xs.group(), subst, no_x, count=1) # only replace first instance in each loop
        no_x = re.sub(r'J', r'X', no_x) # replace J (from unknown translations) with X as these are genuine unknowns
        return no_x

    def match(self, seq, contig, start, end):
        length = len(seq) -6
        pattern = start + r'.' * length + end # regex for matching
        matches = []
        for sequence in contig.sixframe.keys():
            match = re.search(pattern, contig.sixframe[sequence]) # search all protein translations for regex
            if match != None:
                matches.append(match)
        return matches

    # writes hsp information
    def printer (self): 
        logfile.write('Hsp ' + str(self.num) +'\n')
        logfile.write(self.aa_seq + '\n')
        logfile.write('Strand: ' + str(self.strand) + '\n')
        logfile.write('Frame: ' + str(self.frame) + '\n')
        logfile.write('Start: ' + str(self.query_start) + ' End: ' + str(self.query_end) + '\n')

#######################################

# Tile_Path
class Tile_Path(object):
    #Initialise taking values from hsp with highest score and parent contig
    def __init__(self, best_hsp, contig): 
        #Positions are of tiled construct relative to contig (query)
        self.start = best_hsp.query_start
        self.end = best_hsp.query_end
        self.aa_seq = best_hsp.aa_seq 
        self.nt_seq = contig.sequence[self.start-1 : self.end] # nt sequence of tile from position on contig
        self.hsps = deque([best_hsp]) # double-ended queue allows insertion of hsps at either end
        self.strand = best_hsp.strand
        self.frame = best_hsp.frame
        self.added = True # Has anything been added this iteration? Initialise as true to allow loop to start! 
        self.tile = False # Has anything other than the first hsp been added to the tile?

    # adds an hsp to the tile path
    def add_hsp(self, hsp, location, fill=0): 
        if location == 'upstream': # add hsp to 5' end
            self.start = hsp.query_start # adjust nucleotide coordinates relative to contig
            self.hsps.appendleft(hsp) # append hsp to list of those included at correct end
            self.nt_seq = contig.sequence[self.start-1: hsp.query_end+fill] + self.nt_seq
        elif location == 'downstream': # add hsp to 3' end
            self.end = hsp.query_end
            self.hsps.append(hsp)
            self.nt_seq = self.nt_seq + contig.sequence[hsp.query_start-fill-1: self.end]
        else: # something is wrong, exit here
            raise SystemExit ('Error: Cannot add hsp to tile')
        self.aa_seq = translate(self.nt_seq, self.strand)
        hsp.used = True # mark hsp as used
        self.added = True # something has been incorporated into tile
        self.tile = True

    # adds an hsp where it overlaps the existing tile
    def add_overlap(self, hsp, location, distance):
            logfile.write ("Hsp number " + str(hsp.num) + " is " + location + " of the tile and overlapping by " + str(distance+1) + "nucleotides\n")
            if (distance +1) % 3 == 0:
                self.add_hsp(hsp, location, 0-(distance+1))
            else:
                logfile.write("Hsp number " + str(hsp.num) + " not in same frame as tile - correcting\n")
                while (distance +1) % 3 !=0:
                    distance = distance-1
                self.add_hsp(hsp, location, 0-(distance+1))
            logfile.write("Overlapping HSP number " + str(hsp.num) + " added " + location + " of tile\n")
            self.printer()

    # adds an hsp where there is an insertion in the contig between the current hsp and the existing tile
    def add_insert(self, hsp, location, distance): 
        difference = distance -1
        if difference == 1 or difference == 2: # assume incorrect insertion - delete inserted nts
            self.add_hsp(hsp, location)
        elif difference % 3 == 0: # assume correct insertion - include nts and translate
            self.add_hsp(hsp, location, difference)
        else: # insertion with frameshift - attempt to fix
            self.correct_frame(hsp, location, distance) 
        logfile.write("Insertion in contig relative to hit\n")
        logfile.write("Hsp number " + str(hsp.num) + " added " + location + " of tile\n")
        self.printer()

    # adds an hsp where there is a substitution or deletion in the contig between the current hsp and the existing tile
    def add_subst(self, hsp, location, distance):
        difference = distance -1
        if difference == 0 or difference % 3 == 0: # deletion or substitution respectively - include and translate nts if present
            self.add_hsp(hsp, location, difference)
        else: # substitution with frameshift
            self.correct_frame(hsp, location, distance)
        logfile.write("Deletion or substitution in contig relative to hit\n")
        logfile.write("Hsp number " + str(hsp.num) + " added " + location + " of tile\n")
        self.printer()

    # makes a correction between hsps when they are not in the same frame
    def correct_frame(self, hsp, location, distance):
        difference = distance -1
        while difference % 3 != 0: # remove nts until difference is divisible by 3
            difference = difference -1
        self.add_hsp(hsp, location, difference) # add hsp with the cropped nts and translate
        logfile.write("Hsp number " + str(hsp.num) + " not in same frame as tile - correcting\n")

    # write tile information
    def printer(self): 
        logfile.write('Sequence: ' + self.nt_seq + '\n')
        logfile.write('Start: ' + str(self.start) + ' End: ' + str(self.end) + '\n')
        logfile.write('Protein: ' + self.aa_seq + '\n')

#######################################
#Global functions

# translates and sequence in the first frame, forward (default) or reverse
def translate (sequence, strand=1): 
    if strand == -1:
        sequence = revcomp(sequence)
    protein = codons (sequence)
    return protein

# returns a protein string for a sequence in frame 1
def codons (sequence):
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


#returns the reverse complement of a sequence
def revcomp (sequence): 
    sequence = sequence[::-1] #reverse sequence  
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    bases = []
    for s in sequence:
        if s.upper() in complement.keys():
            bases.append(complement[s.upper()])
        else:
            bases.append('N') # can't complement unknown bases so substitute N
    revcomp = ''.join(bases)
    return revcomp


#######################################
# Global Variables 
# Parsing command line options
parser = argparse.ArgumentParser(description="Takes a fasta file of sequences and a BLASTX annotation of that file in xml format.  Attempts to tile Hsps for the highest scoring hit for each sequence, correcting frameshifts in order to improve subsequent annotations.")
# name of fasta file 
parser.add_argument("-f", "--fasta", required=True, help="Fasta file containing sequences")
# name of annotation file (in xml format as code currently stands)
parser.add_argument("-a", "--annotation", required=True, help="Blastx xml file containing annotations for sequences")
# gap limit
parser.add_argument("-g", "--gap_limit", type=int, default=15, help="Cutoff for distance between hsps. If the gap between hsps in nucleotides is greater than this value, the hsp will not be added to the tile.  Default = 15nt")
# print help message if no arguments are given
if len (sys.argv) == 1: 
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

# check that files open
try:
    assembly = open(args.fasta).read().rstrip('\n')
    xml_tree = Blast_Handler(ET.parse(args.annotation))
except IOError:
    raise SystemExit ("Could not open files - please check and try again")
    

#######################################
#Main

# open output files
logfile = open('hsp_tiler_log.txt', 'wt')
outfile = open('hsp_tiler_output.fasta', 'wt')

# split contig file to create a list of contigs
sequences = assembly.split('>')

# create a Contig instance for each contig and append to list
contigs = [] 
for contig in sequences:
    if contig == '':
        continue # ignore empty entry at beginning of list
    else:
        contigs.append(Contig(contig))

# get list of hsps for highest scoring hit for each contig - go to next contig if no hits
for contig in contigs:
    logfile.write('\nNext contig:\n')
    contig.printer()
    xml_tree.get_hsps(contig.name)
    if xml_tree.hsps == None:
        logfile.write('No hits for ' + contig.name + '\n')
        outfile.write('>' + contig.name + '\n' + contig.sequence +'\n')
        continue
    else:
        hsp_list = []
        for each_hsp in xml_tree.hsps: # for each hsp creat hsp instance and add to list
            hsp_list.append(Hsp(each_hsp, contig))
         
        # Initialise tile using highest scoring hsp (first in list) 
        tile = Tile_Path (hsp_list[0], contig) 
        hsp_list[0].used = True # mark that this hsp has been used
        logfile.write('Starting tile...\n')
        hsp_list[0].printer()
        tile.printer()

        # go to next contig if only one hsp in list
        if len(hsp_list) == 1:
            logfile.write('Only one hsp for this hit.  Tile is complete\n')
            outfile.write('>' + contig.name + '\n' + contig.sequence + '\n')
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
                            logfile.write("Hsp number " + str(hsp.num) + "is not on the same strand as the tile - cannot tile\n")
                            hsp.used = True
                            continue # next hsp

                        # ignore hsps that are completely within the tile
                        if tile.start <= hsp.query_start and hsp.query_end <= tile.end: 
                            logfile.write('Hsp ' + str(hsp.num) + ' completely within tile -skipped\n')
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
                            if distance > args.gap_limit: 
                                logfile.write('Too far away to tile\n')
                            # go to next hsp but do not mark as incorporated as may be able to add in a later iteration    
                                continue 

                            # if contig is close enough to attempt tiling                
                            else: 
                                logfile.write("Hsp " + str(hsp.num) + " is " + str(distance -1) + " nucleotides " + location + " of tile\n")
                                assert distance -1 >=0, 'Cannot tile hsps' #sanity check

    # check whether the gap between the hsps relative to the HIT is 0.  If this is true but there is a nucleotide gap relative to the contig, something must have been inserted into the contig.  Statement is unwieldy because test must be made for every combination of hsps upstream and downstream of the tile and on positive and negative strands.
                                if (location == 'upstream' and ((tile.strand == 1 and tile.hsps[0].hit_start - hsp.hit_end -1 == 0) or (tile.strand == -1 and hsp.hit_start - tile.hsps[0].hit_end -1 ==0))) or (location == 'downstream' and ((tile.strand == 1 and hsp.hit_start - tile.hsps[-1].hit_end -1 == 0) or (tile.strand == -1 and tile.hsps[0].hit_start - hsp.hit_end -1 == 0))): 
                                    tile.add_insert(hsp, location, distance)

    #  if the above is not true, check whether the gap between the hsps relative to the hit is > 0.  If this is true and there is no nucleotide gap relative to the contig, something must have been deleted from the contig.  If this is true and there is a nucleotide gap, the nucleotides in that gap must have been substituted in the contig. Statement is unwieldy as above
                                elif (location == 'upstream' and ((tile.strand == 1 and tile.hsps[0].hit_start - hsp.hit_end -1 > 0) or (tile.strand == -1 and hsp.hit_start - tile.hsps[0].hit_end -1 > 0))) or (location == 'downstream' and ((tile.strand == 1 and hsp.hit_start - tile.hsps[-1].hit_end -1 > 0) or (tile.strand == -1 and tile.hsps[0].hit_end - hsp.hit_start -1 > 0))): 
                                    tile.add_subst(hsp, location, distance)
                               
    # as we have already dealt with hsps that overlap the tile, if neither of the above are true, there is a problem somewhere
                                else:
                                    logfile.write('Contig ' + contig.name + ': HSP' + str(hsp.num) + ' could not be tiled')
                                    hsp.used = True
                                    continue

        # add unchanged contig to output if no changes made.  Add tile sequence if changes have been made.
        if tile.tile == False:
            outfile.write('>' + contig.name + '\n' + contig.sequence + '\n')
        elif tile.tile ==True:
            outfile.write('>' + contig.name + '\n' + tile.nt_seq + '\n')

logfile.close()
outfile.close()

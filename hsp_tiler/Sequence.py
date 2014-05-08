# Author: Eli Draizen, Kathryn Crouch
# Date: 06-05-2014
# File: Sequence.py

#Standard libraries
import string

class TranslationMixins:
    """Base class for converting a sequence to protein
    """
    def translate_sequence(self, strand=1, start=0, stop=None):
        """Translates and sequence in the first frame, 
        forward (default) or reverse.

        Parameters:
        ___________
        strand : 1 for forward (default); -1 for reverse

        Return:
        _______
        protein - amino acid sequence
        """
        if strand == 1:
            protein = self.codons(start=start, stop=stop)
        else:
            protein = self.revcomp().codons(start=start, stop=stop)
        
        return protein

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
    def codons(self, start=0, stop=None):
        """Convert nucleotide sequence in frame 1 to a protein sequence

        Parameters:
        ___________
        start : index to start sequence
        stop : index to stop sequence

        Returns:
        ________
        protein : protein string for a sequence in frame 1
        """
        if stop is None:
            stop = len(self)

        aas = []
        for index in xrange(start, stop, 3):
            codon = self.sequence[index:index+3]
            if len(codon) == 3: # 'codons' that aren't 3 bases should only occur at the end of a sequence - can be discarded
                if codon.upper() in TranslationMixins.gencode:
                    aas.append(TranslationMixins.gencode[codon.upper()])
                else:
                    aas.append('J') # substitute J if triplet can't be translated - X not used to prevent unending loops in replace_x
        protein = ''.join(aas)
        return Sequence(name=self.name,
                        sequence=protein)

    complement_table = string.maketrans("ACGT", "TGCA")
    def revcomp(self):
        """returns the reverse complement of a sequence. Adapted from Kevin Karplus'
        BME 205 assignment 1 at UCSC.

        Return:
        ___________
        revcomp : string containing the reverse complement of thenucleotide sequence
        """
        return DNASequence(name="{}_rc".format(self.name),
                           sequence=self.sequence[::-1].translate(TranslationMixins.complement_table)) 

    def sixframe(self):
        """Calculate the six frame translation of a nucleotide sequence and yield
        the frame and preotein sequence.

        Returns:
        ________
        frame : current frame
        protein : protein sequence of the current frame 
        """

        for strand, seq in ((1, self), (-1, self.revcomp())):
            for i in xrange(0, 3):
                protein = seq.codons(start=i)
                frame = (i+1)*strand
                protein.name += "_f{}".format(frame)
                yield frame, protein

class Sequence(object):
    def __init__(self, name="", description="", sequence=""):
        #Initialize an empty sequence
        self.name = name
        self.description = description
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        space = len(self.description)>0 and " " or ""
        tmp = ">{}{}{}\n".format(self.name, space, self.description)
        for i in range(0, len(self.sequence), 60):
            tmp += "{}\n".format(self.sequence[i:i+60])
        return tmp[:-1]

class DNASequence(Sequence, TranslationMixins):
    pass

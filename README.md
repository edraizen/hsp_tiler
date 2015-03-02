# hsp_tiler
A python script to mend frameshift mutions from blast annotation

RNA sequencing (RNA-seq) has become a popular way to study the transcriptome, RNA transcripts present in cell at given time, using Next-Generation Sequencing technologies. RNA-seq is used to quantify gene expression, single nucleotide polymorphism, and gene fusion. Once the RNA has been sequenced, it must be assembled into longer reads, or contigs. This process is simple if there is a reference genome, but becomes more challenging without one, which we call de novo assembly.

The first step in a de novo transcriptome assembly is to annotate the contigs with their functionality. RNA-seq analysis fails if there are any errors in sequences due to sequencing error such as frameshift muta- tions. If not corrected, the errors cause further annotation to fail (i.e. The answer to ‘What is the function of the gene product?’ will be different if frameshifts are present). Here I present HSP-Tiler, a new frameshift correction tool for de novo RNA-seq transcript assemblies, using a combination of BLASTX annotation for similarity-based error correct and a hidden Markov Model (HMM) for Ab initio error correction.

The program requires two inputs, a FASTA file and a BLASTX annotation file, but adding other parameters described below may produce better results:

• -f, –fasta Path to FASTA file containing sequences which are thought to have frameshift mutations.

• -a, –annotation Path to annotation file, which is by default a BLASTX output in BLAST-XML format. Each contig query will have a hits containing HSPs, or high-scoring segment pairs. HSP-Tiler does not run BLAST itself, so the user must run BLASTX before using HSP-Tiler. BLASTX against the non-redundant (nr) database from NCBI leads to the best results.

• –format Format of the annotation file. Can use output from BLASTX (xml (default), or tabs with specified format), BLAT PSLX, and RapSearch2.

• -g, –gap_limit Tell the program to only include HSPSs where the gap between the current HSP and the tile are less than the gap cutoff, otherwise they are ignored.

• -e, –evalue_cutoff Tell the program to only include an HSP if the e-value is less than or equal to cutoff, else they are ignored.

• –allHits Use all of the hits from BLASTX, instead of just the 1st one. Default is false.

• –filter Tell the program to only use hits from genus or species who are not part of the filter species or genus. Example: a filter of ‘Caenorhabditis’ excludes all hits from the genus Caenorhabditis, and ‘Caenorhabditis elegans’ will only exclude the species Caenorhabditis elegans.

• –filterType Regex to retrieve taxon info from hit. Use 0 for NCBI, or 1 for uniprot, or the full regex query

• –codon_usage Path to file, URL, or ID of file containing codon usage.

• –compute_codon Compute the codon usage table for given FASTA sequences.

• -o, –outfile File to save corrected sequences.

• -l, –logfile File to save log.

• –ignoreNoHit Ignore sequence if there are no BLAST hits

• –noHitORF Output the longest open reading frame for the sequences with no BLAST hits

• –noHitFrame1 Output the first reading frame for the sequences with no BLAST hits • -p, –protein Output protein sequence. Default is false.

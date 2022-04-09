 # SPINE: Saturated Programmable INsertion Engineering
### Protein domain insertion via programmed oligo libraries
Python script for generating oligo libraries and PCR primers for programmed domain insertion

# Python Dependencies
biopython <br />
numpy

# Data Formats
## Input:
-> Targeted genes must be in fasta format and include a minimum of 30 bases surrounding gene.
-> Entire plasmid sequence with all the essential genetic context is advised to search for nonspecific amplification <br />
-> Multiple sequences can be input for deep mutation scanning in one file.

## Output:
-> Final output is fasta format. One file for oligo pools and one file for PCR primers
-> Different mutation type function can be run separately. The barcode might be overlapped if the functions are run separately. "allmut" function wrapped up everything and generate oligo libraries based on the same fragmentation strategy of the GOI. The primers for both of the oligo libraries or plasmid vector within the same fragmentation part are generic.

# Notes:
Gene primers should be same for the same sequence, because no matter what the mutation types we introduced here the fragmentation strategy should always be the same based on the size of the GOI

# Position arguments
Gene start is defined as base number of first base in first codon and gene end is defined as base number of last base in last codon which should be the stop codon.
(Program will subtract 1 from gene start for python numbering)


To define gene position within given fasta file, add start:# end:# to fasta description. (Otherwise use command line prompts) <br />
'>geneA start:11 end:40'


# Pairs of files for running
-> Final version for running DMS: 'main_spine.py' with functional module 'SPINE_v7.py' 
	SPINE_v7 contains all of the single amino acid mutation function and also the integrated version for bulk deep mutation scanning
		
		DMS: single aa replacement
		
		S_INS/S_DEL: single aa insertion / deletion
		
		allmut: combination of all single aa mutation types
		
		allmut_noCys: exclude introduction of Cysteine from the all mutation types
	
	SPINE_v7 also allows for multiple input of GOI sequences in one input file
	
	SPINE_v7 carries counter to generate the diversity of final oligo libraries



# Running Test
Domain Insertion Scanning:
Note: Domain Insertion scanning is directly from the original script. I did not make any changes.

	python3.8 main_spine.py -wDir tests -geneFile combined_fasta.fa -oligoLen 230 -mutationType DIS

Deep Mutational Scanning:

	python3.8 main_spine.py -wDir test_bulk_mut -geneFile Multiple_inputfile.fa -oligoLen 230 -mutationType allmut -usage ecoli

	python3.8 main_spine.py -wDir test_bulk_mut -geneFile Kir.fa -oligoLen 230 -mutationType DMS -usage ecoli



# Usage
```
optional arguments:
-h, --help                 show this help message and exit
-wDir WDIR                 Working directory for fasta files and output folder
-geneFile GENEFILE         Input all gene sequences including backbone in a fasta
                           format. Place all in one fasta file. Name description
                           can include start and end points (>gene1 start:1
                           end:2)
-handle HANDLE             Genetic handle for domain insertion.  This is important
                           for defining the linker. Currently uses BsaI (4 base
                           overhang), but this can be swapped for SapI (3 base
                           overhang).
-matchSequences            Find similar sequences between genes to avoid printing
                           the same oligos multiple times. Default: No matching
-oligoLen OLIGOLEN         Synthesized oligo length
-fragmentLen FRAGMENTLEN   Maximum length of gene fragment
-overlap OVERLAP           Enter number of bases to extend each fragment for
                           overlap. This could help with insertion coverage close to
                           fragment boundary. Overlap does not add additional
                           insertions and thus no additional oligos.
-mutationType              Run deep insertion scan "DIS", deep mutation scan for single amino acid replacement "DMS", single amino acid deletion "S_DEL", 
			  single amino acid insertion "S_INS" or combine all single amino acid mutation together(replacement, deletion, insertion) "allmut".
			  Or exclude any introduction of 'Cys' into the mutation library "allmut_noCys" 
-usage USAGE               Default is "human". Or select "ecoli"
```

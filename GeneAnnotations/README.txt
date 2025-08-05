This folder contains the scripts used for gene annotation analyses in the associated publication https://doi.org/10.1093/gbe/evac090

Scripts were used as follows:
Wilcox-IndexRef(FalcoRef8).bash simply indexed reference genomes for these analyses
Wilcox-Fasta2Phylip(FalcoRef8).pl was used to convert fasta to relaxed phylip format. This script was not written by me and was taken from https://github.com/npchar/Phylogenomic/blob/master/fasta2relaxedPhylip.pl
Wilcox-AlignCDS(FalcoRef8).bash performed codon aware alignments of CDS regions with PRANK and creates input files for PAML
Wilcox-ReAnnotateAlignCDS(FalcoRef8).bash reannotated Phylip alignments for use in PAML

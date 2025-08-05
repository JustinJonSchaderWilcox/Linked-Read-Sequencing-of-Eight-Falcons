#Converting output of RAxML tree of falcon genome variables sites to Newick format
perl -pe 's#(:0\.[0-9]+)\[([0-9]+)+\]#$2$1#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/RAxML_bipartitionsBranchLabels.Ref8T3 > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ref8Gamma_RAxML_bipartitionsBranchLabels.tre
cp /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/RAxML_bipartitionsBranchLabels.Ref8T3 /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ref8Gamma_RAxML_bipartitionsBranchLabels.nwk

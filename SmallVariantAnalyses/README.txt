This folder contains the scripts used for analyses using small variants in the associated publication https://doi.org/10.1093/gbe/evac090

The scripts were used as follows:
Wilcox-LongRangerAlignFalco(FalcoRef8).bash aligned raw 10X linked reads to the haploid assemblies of the large falcons
Wilcox-LongRangerAlignKestrel(FalcoRef8).bash aligned raw 10X linked reads to the haploid assembly of the common kestrel
Wilcox-MUMmerKestrelAlignAllGenomes(FalcoRef8).bash aligned all assemblies to the kestrel assembly
Wilcox-MUMmerSelfAlignAllGenomes(FalcoRef8).bash aligned all phased haplotype assemblies (pseudohaplotypes) to their respective complement
Wilcox-MUMmer2VCFkestrelAlign(FalcoRef8).bash converted MUMmer outputs from alignments to the kestrel assembly to VCFs
Wilcox-MUMmer2VCFselfAlign(FalcoRef8).bash converted MUMmer outputs from alignments of pseudohaplotypes to VCFs
Wilcox-PhasedKestrelCorrectHetIntersects(FalcoRef8).bash called heterozygous sites within kestrel MUMmer pseudohaplotype alignments
Wilcox-IntersectSelfAlignAllvLongRangerGenomesVCFs(FalcoRef8).bash called heterozygous sites between aligned pseudohaplotypes
Wilcox-VariantSelfAlignCompile(FalcoRef8).bash compiled small variant counts across sources for each individual
Wilcox-KestrelAlignUltimateCompilation(FalcoRef8).bash compiled variants counts and types across methods from alignments to the kestrel assembly
Wilcox-CompileBaseCompositionConservedSites(FalcoRef8).bash compiled base composition at conserved/invariant sites
Wilcox-CompileLengthsUniqIndels(FalcoRef8).bash compiled length information for indels that were private to individuals
Wilcox-KestrelCongnAlignVSconMA(FalcoRef8).bash created consensus genomes of falcon based on alignments to kestrel genome to produce alignments of only variable sites for phylogenetic analysis
Wilcox-AllFalcoRAxML(FalcoRef8).bash ran RAxML on all falcon genomes to create a maximum-likelihood phylogenetic tree
Wilcox-RAxML2Tre(FalcoRef8).bash converted RAxML tree formats (Newick) to '.tre' (Nexus) formats
Wilcox-AllFalcoUltrametricPhylo(FalcoRef8).bash created an ultrametric time-scaled tree
Wilcox-PhyloCladeSharedVar(FalcoRef8).bash compiled small variants types by tree branches
Wilcox-AllFalcoPSMC(FalcoRef8).bash ran PSMC on falcon assemblies
Wilcox-AllFalcoBootStrapPSMC(FalcoRef8).bash ran bootstrapped PSMC analysis on falcon assemblies

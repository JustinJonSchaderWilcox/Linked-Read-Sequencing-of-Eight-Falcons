This folder contains the scripts used for repeat and structurla variant analyses in the associated publication https://doi.org/10.1093/gbe/evac090

The scripts were used as follows:
Wilcox-AllFalcoPhasedRepeatModeler2(FalcoRef8).bash ran RepeatModeler2 on diploid falcon assemblies
Wilcox-AllFalcoMergeRepeatModeler2(FalcoRef8).bash merged and clustered RepeatModeler2 outputs for use in Repeatmasker
Wilcox-AllFalcoSizeFilterMergedRepeatModeler2(FalcoRef8).bash filtered merged and clustered RepeatModeler2 outputs to at least 40 copies across all falcons
Wilcox-AllFalcoRepeatModeler2ConsensusLibrary(FalcoRef8).bash built a consensus library of merged and clustered RepeatModeler2 libraries for use with RepeatMasker
Wilcox-AllFalconRepMaskRMBLAST(FalcoRef8).bash ran RepeatMasker with rmblast using merged and clustered consensus libraries from RepeatModeler2
Wilcox-RepeatCompareOtherGenomesDownload(FalcoRef8).bash downloaded from other birds for comparison of repeat landscapes
Wilcox-RepmaskOtherGenome(FalcoRef8).bash ran RepeatMasker on other avian genomes using DFAM and HMMER
Wilcox-RepMod2OtherGenome(FalcoRef8).bash ran RepeatModeler2 on other avian genomes
Wilcox-RepmaskRepmod2OtherGenome(FalcoRef8).bash ran RepeatMasker using de novo RepeatModeler2 libraries on other avain genomes
Wilcox-RepMaskGraphOtherGenome(FalcoRef8).bash creates divergence curves for other avian genomes
Wilcox-RepSVannotate(FalcoRef8).bash annotated large structural variants in genomes based on MUMmer alignments (see 'SmallVariantAnalyses') and checks for intersections with repeats
Wilcox-SummarizeSVs(FalcoRef8).bash compiled and summarized large structural variant outputs
Wilcox-PairRepSVannotate(FalcoRef8).bash assessed polymorphic variation in structural variants using pairwise comparisons between falcon psuedohaplotypes
Wilcox-RepeatORF(FalcoRef8).bash annotated ORFs within annotated repeats
Wilcox-RepeatActiveORF(FalcoRef8).bash searched for potentially active ORFs within annotated repeats of falcons

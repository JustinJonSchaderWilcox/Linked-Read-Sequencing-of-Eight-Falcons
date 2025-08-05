#(re)Compiling the base composition of conserved sites
##Conserved Sites Summary
srun --pty -n 1 /bin/bash
A=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/vReference/0000.vcf | awk 'BEGIN{Count=0}$4=="A"{Count=Count+1}END{print Count}'`
T=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/vReference/0000.vcf | awk 'BEGIN{Count=0}$4=="T"{Count=Count+1}END{print Count}'`
G=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/vReference/0000.vcf | awk 'BEGIN{Count=0}$4=="G"{Count=Count+1}END{print Count}'`
C=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/vReference/0000.vcf | awk 'BEGIN{Count=0}$4=="C"{Count=Count+1}END{print Count}'`
CombinedAT=`expr $A + $T`
CombinedGC=`expr $G + $C`

rm /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/HaploidConservedSitesBaseComp.txt
echo -e Conserved_Bases'\t'A'\t'T'\t'G'\t'C'\t'CombinedAT'\t'CombinedGC > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/HaploidConservedSitesBaseComp.txt
echo -e Haploid_Count'\t'$A'\t'$T'\t'$G'\t'$C'\t'$CombinedAT'\t'$CombinedGC >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/HaploidConservedSitesBaseComp.txt

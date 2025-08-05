#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --mem=118GB
#SBATCH -p serial
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=28
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out


rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric

echo "#NEXUS" > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.ultrametric.rs8_ctl
echo "begin trees;" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.ultrametric.rs8_ctl
echo "tree Ref8 =" $(cat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ref8Gamma_RAxML_bipartitionsBranchLabels.tre) >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.ultrametric.rs8_ctl
echo "End;" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.ultrametric.rs8_ctl
echo "begin rates;" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.ultrametric.rs8_ctl
echo "blformat lengths=persite nsites=28479719 ultrametric=no round=yes;" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.ultrametric.rs8_ctl
echo "collapse;" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.ultrametric.rs8_ctl
echo "mrca Kestrel A_genome Kest2_geno;" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.ultrametric.rs8_ctl
echo "fixage taxon=Kestrel age=7246;" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.ultrametric.rs8_ctl
echo "divtime method=pl algorithm=tn cvStart=0 cvInc=0.5 cvNum=8 crossv=yes;" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.ultrametric.rs8_ctl
echo "describe plot=chronogram;" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.ultrametric.rs8_ctl
echo "describe plot=tree_description;" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.ultrametric.rs8_ctl
echo "end;" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.ultrametric.rs8_ctl
/project/jw5478/Programs/r8s1.81/src/r8s -b -f /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.ultrametric.rs8_ctl > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.r8s
tail -1 /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.r8s > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.tre
sed -i 's#A_genome#Lanner#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.tre
sed -i 's#B_genome#Barbay#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.tre
sed -i 's#C_genome#Peregrine#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.tre
sed -i 's#D_genome#Saker#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.tre
sed -i 's#E_genome#Gyr-1#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.tre
sed -i 's#F_genome#Gyr-2#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.tre
sed -i 's#FpBSx77_ge#Bl. Shaheen#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.tre
sed -i 's#Kest2_geno#Kestrel#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.tre
sed -i 's#Kestrel;#;#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.tre

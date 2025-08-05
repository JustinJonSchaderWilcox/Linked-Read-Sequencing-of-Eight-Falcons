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

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML
cd /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML
cp /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Alignments/MA_Kest100KBhap.VarSite_dip_con.vs.phylip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/MA_Kest100KBhap.VarSite_dip_con.vs.phylip
/home/jw5478/Executables/Programs/standard-RAxML/raxmlHPC-PTHREADS -s /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/MA_Kest100KBhap.VarSite_dip_con.vs.phylip -p 31415927 -o Kest2_geno -n Ref8T1 -m GTRGAMMA -T 28
/home/jw5478/Executables/Programs/standard-RAxML/raxmlHPC-PTHREADS -s /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/MA_Kest100KBhap.VarSite_dip_con.vs.phylip -p 31415927 -o Kest2_geno -n Ref8T2 -m GTRGAMMA -b 271828 -# 1000 -T 28
/home/jw5478/Executables/Programs/standard-RAxML/raxmlHPC-PTHREADS -s /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/MA_Kest100KBhap.VarSite_dip_con.vs.phylip -f b -t RAxML_bestTree.Ref8T1 -z RAxML_bootstrap.Ref8T2 -p 31415927 -o Kest2_geno -n Ref8T3 -m GTRGAMMA

#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --mem=116GB
#SBATCH -p serial
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=28
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

s=$1
MUMmer=/home/jw5478/Executables/Programs/MUMmer3.23/

cd /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/MUMmerOutputs/ 

$MUMmer/nucmer --maxgap=500 --mincluster=100 /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta  /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_16KBph2.diploid.fasta --prefix="$s"-v-"$s"100KBhap
$MUMmer/show-coords -c -d -l -r -T /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/MUMmerOutputs/"$s"-v-"$s"100KBhap.delta > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/MUMmerOutputs/"$s"-v-"$s"100KBhap.coords
$MUMmer/show-snps -r -T /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/MUMmerOutputs/"$s"-v-"$s"100KBhap.delta > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/MUMmerOutputs/"$s"-v-"$s"100KBhap.snps

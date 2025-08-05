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


cd /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/

$MUMmer/nucmer --maxgap=500 --mincluster=100 /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Falco_rusticolus.chrom.fna /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta --prefix="$s"-v-GyrChrom
$MUMmer/show-coords -c -d -l -r -T /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-GyrChrom.delta > /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-GyrChrom.coords
$MUMmer/show-snps -r -T /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-GyrChrom.delta > /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-GyrChrom.snps

$MUMmer/nucmer --maxgap=500 --mincluster=100 /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Falco_peregrinus.pseudo_chrom.fna /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta --prefix="$s"-v-pseudoPeregrineChrom
$MUMmer/show-coords -c -d -l -r -T /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-pseudoPeregrineChrom.delta > /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-pseudoPeregrineChrom.coords
$MUMmer/show-snps -r -T /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-pseudoPeregrineChrom.delta > /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-pseudoPeregrineChrom.snps

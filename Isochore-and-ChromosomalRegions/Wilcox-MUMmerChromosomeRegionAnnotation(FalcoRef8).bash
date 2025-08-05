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

$MUMmer/nucmer --maxgap=500 --mincluster=100 /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Thrush.chrom.fna /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta --prefix="$s"-v-ThrushChrom
$MUMmer/show-coords -c -d -l -r -T /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-ThrushChrom.delta > /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-ThrushChrom.coords
$MUMmer/show-snps -r -T /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-ThrushChrom.delta > /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-ThrushChrom.snps

$MUMmer/nucmer --maxgap=500 --mincluster=100 /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Kakapo.chrom.fna /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta --prefix="$s"-v-KakapoChrom
$MUMmer/show-coords -c -d -l -r -T /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-KakapoChrom.delta > /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-KakapoChrom.coords
$MUMmer/show-snps -r -T /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-KakapoChrom.delta > /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-KakapoChrom.snps

$MUMmer/nucmer --maxgap=500 --mincluster=100 /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Seriema.chrom.fna /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta --prefix="$s"-v-SeriemaChrom
$MUMmer/show-coords -c -d -l -r -T /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-SeriemaChrom.delta > /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-SeriemaChrom.coords
$MUMmer/show-snps -r -T /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-SeriemaChrom.delta > /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-SeriemaChrom.snps

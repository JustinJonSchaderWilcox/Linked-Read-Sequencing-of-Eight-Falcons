#!/bin/bash
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=28
#SBATCH --mem=250GB
# Output and error files
#SBATCH -o job.%J.out
#SBATCH --time=168:00:00
#SBATCH --partition=c18_25

#Modules
module purge
module load perl/intel/5.24.0

#Input
s=$1
RepMod="/home/jw5478/Programs/RepeatModeler-2.0.1"

#Script
##Database

rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepModDB/"$s"
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepModDB/"$s"

cd /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepModDB/"$s"

#Build Database
$RepMod/BuildDatabase -engine ncbi -name "$s"_FalcoRefGenomeRepModDB /scratch/jw5478/FalconProject/Projects/RefGenomes/Genomes/SuperNova10X/16KBph2Diploid/"$s"_16KBph2.diploid.fasta

##RepeatModeler
$RepMod/RepeatModeler -pa 7 -database "$s"_FalcoRefGenomeRepModDB -engine ncbi -genomeSampleSizeMax 243000000 -LTRStruct

#Cleaning up intermediary RepeatModeler files
rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepModDB/"$s"/RM_*

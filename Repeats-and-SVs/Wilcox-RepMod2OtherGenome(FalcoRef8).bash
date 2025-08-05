#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180GB
#SBATCH --partition=cs
#SBATCH -o job.%J.out
#SBATCH --time=168:00:00

#Modules
module purge

#Varialbes
o=$1

GENOME=`echo $o | perl -pe 's#Thrush#GCA_009819885.1_bCatUst1.pri_genomic.fna#' | perl -pe 's#Kakapo#GCF_004027225.2_bStrHab1.2.pri_genomic.fna#' | perl -pe 's#Seriema#GCA_009819825.1_bCarCri1.pri_genomic.fna#'`

#Input
RepMod=/home/jw5478/Programs/RepeatModeler-2.0.1
rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/"$o"
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/"$o"

#Script
cd /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/"$o"
#Build Database
$RepMod/BuildDatabase -engine ncbi -name "$o"_OtherGenomeRepModDB /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/Genomes/"$o"/"$GENOME"
##RepeatModeler
$RepMod/RepeatModeler -pa 12 -database "$o"_OtherGenomeRepModDB -engine ncbi -genomeSampleSizeMax 123000000 -LTRStruct

#Cleaning up intermediary RepeatModeler files
rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/"$o"/RM_*

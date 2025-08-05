#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180GB
#SBATCH -o job.%J.out
#SBATCH --time=168:00:00

#Modules
module purge
module load busco/5.0.0
module load augustus/intel/3.4.0

#Variable
s=$1

rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Genes/BUSCO/"$s"
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Genes/BUSCO/"$s"
cp  /scratch/jw5478/FalconProject/Projects/RefGenomes/Genomes/SuperNova10X/16KBph2Haploid/"$s"_16KBph2.haploid.fasta /scratch/jw5478/FalconProject/Projects/RefGenomes/Genes/BUSCO/"$s"/"$s"_16KBph2.haploid.fasta

cp -r /share/apps/augustus/3.4.0/intel/config /scratch/jw5478/FalconProject/Projects/RefGenomes/Genes/BUSCO/"$s"/config
export AUGUSTUS_CONFIG_PATH="/scratch/jw5478/FalconProject/Projects/RefGenomes/Genes/BUSCO/"$s"/config"

cd /scratch/jw5478/FalconProject/Projects/RefGenomes/Genes/BUSCO/"$s"
busco -i ./"$s"_16KBph2.haploid.fasta -o "$s"_BUSCO_Genome -l aves_odb10 --augustus -m genome -c 48 --long --limit 5 --augustus_species chicken

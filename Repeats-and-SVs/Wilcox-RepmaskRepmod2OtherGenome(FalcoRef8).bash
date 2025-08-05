#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180GB
#SBATCH --partition=cs
#SBATCH -o job.%J.out
#SBATCH --time=168:00:00

module purge
module load perl/intel/5.32.0


o=$1
GENOME=`echo $o | perl -pe 's#Thrush#GCA_009819885.1_bCatUst1.pri_genomic.fna#' | perl -pe 's#Kakapo#GCF_004027225.2_bStrHab1.2.pri_genomic.fna#' | perl -pe 's#Seriema#GCA_009819825.1_bCarCri1.pri_genomic.fna#'`


rm /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/"$o"/"$o"_OtherGenomeRepModDB-families.size-label.fa
rm /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/"$o"/"$o"_OtherGenomeRepModDB-families.filt-size20.fa
rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/"$o"
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/"$o"

perl -pe 's#\s\(.*Final Multiple Alignment (Size = [0-9]+).*#;$1#' /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/"$o"/"$o"_OtherGenomeRepModDB-families.fa | perl -pe 's#Size = #size=#' > /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/"$o"/"$o"_OtherGenomeRepModDB-families.size-label.fa
/home/jw5478/Programs/USEARCH/usearch11.0.667_i86linux32 -sortbysize /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/"$o"/"$o"_OtherGenomeRepModDB-families.size-label.fa -fastaout /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/"$o"/"$o"_OtherGenomeRepModDB-families.min-size20.fa -minsize 20

/home/jw5478/Programs/RepeatMasker/RepeatMasker -e rmblast -a -lib /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/"$o"/"$o"_OtherGenomeRepModDB-families.min-size20.fa -dir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/"$o" -gccalc -poly -s -pa 12 /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/Genomes/"$o"/"$GENOME"

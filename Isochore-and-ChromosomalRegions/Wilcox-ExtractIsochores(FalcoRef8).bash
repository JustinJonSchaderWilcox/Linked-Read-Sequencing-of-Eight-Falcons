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

module purge
module load perl
module load gencore/1
module load gencore_variant_detection/1.0
module load gencore_annotation/1.0

mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Support
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Beds
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Seqs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Variants
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Variants/SNVs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Summary
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Summary/Composition

for s in A D B FpBSx77 C E F Kestrel
do
awk 'BEGIN{OFS=FS="\t"}$2>=100000{print $1, $2}' /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta.fai | sort -k1,1V -k2,2n > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Support/"$s"_16KBph2.100KBhaploid.genome
bedtools makewindows -g /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Support/"$s"_16KBph2.100KBhaploid.genome -w 100000 | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2+1 OFS $3}' | awk '{Size=$3-$2}Size==99999{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Beds/"$s"_16KBph2.haploid.win100KB.bed &
done
wait

#Function will extract sequences from these; paths need to be changed to extract from files for this analysis
function extract_regions {
WINDOW_POSITIONS=`cat /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Beds/"$s"_16KBph2.haploid.win100KB.bed | perl -pe 's#\b\s+\b#__#g'`
rm /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Seqs/"$s"_16KBph2.haploid.win100KB.fasta
for p in $WINDOW_POSITIONS
do
REF_SCAFFOLD=`echo $p | perl -pe 's#__#\t#g' | awk '{print $1}'`
REF_POS1=`echo $p | perl -pe 's#__#\t#g' | awk '{print $2}'`
REF_POS2=`echo $p | perl -pe 's#__#\t#g' | awk '{print $3}'`
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta $REF_SCAFFOLD:$REF_POS1-$REF_POS2 >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Seqs/"$s"_16KBph2.haploid.win100KB.fasta
done
}

for s in A D B FpBSx77 C E F Kestrel
do
extract_regions &
done
wait

#Pull out base composition of each Isochore
Seqtk="/home/jw5478/Executables/Programs/seqtk/seqtk"
for s in A D B FpBSx77 C E F Kestrel
do
echo 'Win\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Summary/Composition/"$s"_16KBph2.haploid.win100KB.comp
$Seqtk comp /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Seqs/"$s"_16KBph2.haploid.win100KB.fasta >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Summary/Composition/"$s"_16KBph2.haploid.win100KB.comp &
done
wait

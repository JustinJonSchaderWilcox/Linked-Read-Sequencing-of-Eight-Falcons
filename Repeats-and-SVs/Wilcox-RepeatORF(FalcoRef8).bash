#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --mem=60GB
#SBATCH -p serial
# Set number of nodes to run
# Set number of tasks to run
#SBATCH --ntasks=14
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

module purge
module load perl
module load gencore/1
module load gencore_variant_detection/1.0
module load jdk

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active/Genomes
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active/AnnotateBeds
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active/AnnotateSeqs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active/ORFs

function repeat_seqs {
REPEAT_POSITIONS=`cat /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active/AnnotateBeds/"$s"_16KBph2.diploid.fasta.all_rep.bed | perl -pe 's#\b\s+\b#__#g'`
for p in $REPEAT_POSITIONS
do
REF_SCAFFOLD=`echo $p | perl -pe 's#__#\t#g' | awk '{print $1}'`
REF_POS1=`echo $p | perl -pe 's#__#\t#g' | awk '{print $2}'`
REF_POS2=`echo $p | perl -pe 's#__#\t#g' | awk '{print $3}'`
NAME=`echo $p | perl -pe 's#__#\t#g' | awk '{print $4}'`
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active/Genomes/"$s"_16KBph2.diploid.fasta $REF_SCAFFOLD:$REF_POS1-$REF_POS2 | perl -pe "s%(>.*)%\1#$NAME%" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active/AnnotateSeqs/"$s"_16KBph2.diploid.fasta.all_rep.fasta
done
}

for s in A B C D E F FpBSx77 Kestrel
do
cp /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_16KBph2.diploid.fasta /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active/Genomes/"$s"_16KBph2.diploid.fasta &
done
wait

for s in A B C D E F FpBSx77 Kestrel
do
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active/Genomes/"$s"_16KBph2.diploid.fasta &
done
wait

for s in A B C D E F FpBSx77 Kestrel
do
grep -P "[_0-9+C]+" /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Annotations/"$s"_16KBph2.diploid.fasta.out | awk 'BEGIN{OFS="\t"}{print $5 OFS $6 OFS $7 OFS $11}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active/AnnotateBeds/"$s"_16KBph2.diploid.fasta.all_rep.bed &
done
wait

for s in A B C D E F FpBSx77 Kestrel
do
repeat_seqs &
done
wait

for s in A B C D E F FpBSx77 Kestrel
do
/home/jw5478/Executables/Programs/EMBOSS-6.6.0/bin/getorf -sequence /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active/AnnotateSeqs/"$s"_16KBph2.diploid.fasta.all_rep.fasta -find 3 -outseq /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active/ORFs/"$s"_16KBph2.diploid.fasta.all_rep.ORF.fasta &
done

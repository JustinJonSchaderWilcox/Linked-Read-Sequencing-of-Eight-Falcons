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
module load gencore
module load gencore_annotation/1.0
#Variables
FalconidaeMito="/scratch/jw5478/FalconProject/Projects/RefGenome8/Mitochondria/NUMT"
HMMER="/home/jw5478/Executables/Programs/HMMER/bin"
Mafft="/home/jw5478/Executables/Programs/Mafft"
Emboss="/home/jw5478/Executables/Programs/EMBOSS-6.6.0/bin"
DATE=`ls $FalconidaeMito/Downloads/Sequences/FalconidaeMitoGenome_*.fasta | perl -pe 's#.*FalconidaeMitoGenome_(.*)\.fasta#\1#'`

#Index Downloaded Genomes
samtools faidx $FalconidaeMito/Downloads/Sequences/FalconidaeMitoGenome_"$DATE".fasta
#Alignment Directory
#Align for Consensus
$Mafft/mafft-linsi --nuc --adjustdirectionaccurately --reorder --thread 8 --maxiterate 1000 $FalconidaeMito/Downloads/Sequences/Compiled/AvesSampledFalconidaeMitoGenome_"$DATE".fasta > $FalconidaeMito/Downloads/Sequences/Alignments/AvesSampledFalconidaeMitoGenome_"$DATE".linsi.fasta

#HMMER Time
$HMMER/hmmbuild $FalconidaeMito/Downloads/Sequences/Alignments/AvesSampledFalconidaeMitoGenome_"$DATE".linsi.hmm  $FalconidaeMito/Downloads/Sequences/Alignments/AvesSampledFalconidaeMitoGenome_"$DATE".linsi.fasta
$HMMER/nhmmer --cpu 12 -A $FalconidaeMito/Downloads/Sequences/Alignments/AvesSampledFalconidaeMitoGenome_"$DATE".stk -E 0.000001 --tblout $FalconidaeMito/Downloads/Sequences/Alignments/Records/AvesSampledFalconidaeMitoGenome_"$DATE".orient.txt --incE 0.000001   $FalconidaeMito/Downloads/Sequences/Alignments/AvesSampledFalconidaeMitoGenome_"$DATE".linsi.hmm $FalconidaeMito/Downloads/Sequences/Compiled/AvesSampledFalconidaeMitoGenome_"$DATE".fasta

#Orient Sequences for Alignment
ACCESSIONS=`awk '!/^#/{print $1}' $FalconidaeMito/Downloads/Sequences/Alignments/Records/AvesSampledFalconidaeMitoGenome_"$DATE".orient.txt | sort | uniq`
rm $FalconidaeMito/Downloads/Sequences/AvesSampledFalconidaeMitoGenome_"$DATE".orient.fasta
for a in $ACCESSIONS
do
DESCRIPTION=`awk -v accession=$a 'BEGIN{OFS="\t"}$1==accession{$1=$2=$3=$4=$5=$6=$7=$8=$9=$10=$11=$12=$13=$14=$15="";print $0; exit}' $FalconidaeMito/Downloads/Sequences/Alignments/Records/AvesSampledFalconidaeMitoGenome_"$DATE".orient.txt | perl -pe 's#\s+# #g'`
ORIENT=`awk -v accession=$a 'BEGIN{OFS="\t"}$1==accession && $5>=$7{start=$5-$7; exit}$1==accession && $7>=$5{start=$7-$5; exit}END{print start}' $FalconidaeMito/Downloads/Sequences/Alignments/Records/AvesSampledFalconidaeMitoGenome_"$DATE".orient.txt`
START=`expr $ORIENT + 1`
samtools faidx $FalconidaeMito/Downloads/Sequences/Compiled/AvesSampledFalconidaeMitoGenome_"$DATE".fasta $a:$START > $FalconidaeMito/Downloads/Sequences/AvesSampledFalconidaeMitoGenome_"$DATE".orient_temp.fasta
samtools faidx $FalconidaeMito/Downloads/Sequences/Compiled/AvesSampledFalconidaeMitoGenome_"$DATE".fasta $a:1-$ORIENT >> $FalconidaeMito/Downloads/Sequences/AvesSampledFalconidaeMitoGenome_"$DATE".orient_temp.fasta
awk -v accession=$a -v description="$DESCRIPTION" 'NR==1{print ">"accession" "description; next}!/>/{printf $0}END{printf "\n"}' $FalconidaeMito/Downloads/Sequences/AvesSampledFalconidaeMitoGenome_"$DATE".orient_temp.fasta >> $FalconidaeMito/Downloads/Sequences/AvesSampledFalconidaeMitoGenome_"$DATE".orient.fasta
done
rm $FalconidaeMito/Downloads/Sequences/AvesSampledFalconidaeMitoGenome_"$DATE".orient_temp.fasta

#Align Orient Circular Sequences and Produce Consensus Sequence and then Concatenate to find Alignments that bridge break in Circle
$Mafft/mafft-linsi --nuc --adjustdirectionaccurately --reorder --thread 8 --maxiterate 1000 $FalconidaeMito/Downloads/Sequences/AvesSampledFalconidaeMitoGenome_"$DATE".orient.fasta | awk 'NR==1{print $0; next}/^>/{print "\n"$0}!/^>/{printf $0}END{printf "\n"}' | awk '/^>/{print $0}!/^>/{print$0$0}' > $FalconidaeMito/Downloads/Sequences/AvesSampledFalconidaeMitoGenome_"$DATE".orient.dup.linsi.fasta
$HMMER/hmmbuild $FalconidaeMito/Downloads/Sequences/Alignments/AvesSampledFalconidaeMitoGenome_"$DATE".orient.dup.linsi.hmm $FalconidaeMito/Downloads/Sequences/AvesSampledFalconidaeMitoGenome_"$DATE".orient.dup.linsi.fasta
for s in A B C D E F FpBSx77 Kestrel
do
$HMMER/nhmmer --cpu 12 -A $FalconidaeMito/Downloads/Sequences/Alignments/"$s"_16KBph2.diploid.numts.dup.stk -E 0.000001 --tblout $FalconidaeMito/Downloads/Sequences/Alignments/Records/"$s"_numtsFalconidaeMitoGenome_"$DATE".orient.txt --incE 0.000001 $FalconidaeMito/Downloads/Sequences/Alignments/AvesSampledFalconidaeMitoGenome_"$DATE".orient.dup.linsi.hmm /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_16KBph2.diploid.fasta
done
for s in A B C D E F FpBSx77 Kestrel
do
$HMMER/esl-reformat fasta $FalconidaeMito/Downloads/Sequences/Alignments/"$s"_16KBph2.diploid.numts.dup.stk | awk 'NR==1{print $0; next}/^>/{print "\n"$0}!/^>/{printf $0}END{printf "\n"}' > $FalconidaeMito/Downloads/Sequences/Alignments/"$s"_16KBph2.diploid.numts.dup.fasta
for n in $(grep ">" $FalconidaeMito/Downloads/Sequences/Alignments/"$s"_16KBph2.diploid.numts.dup.fasta | perl -pe 's#(>[0-9]+_[1-2]+/[0-9]+-[0-9]+).*#$1#' | sort -V | uniq)
do
awk -v name="$n" '$0~name{numt="YES"; next}numt=="YES"{print name"\n"$0; exit}' $FalconidaeMito/Downloads/Sequences/Alignments/"$s"_16KBph2.diploid.numts.dup.fasta >> $FalconidaeMito/Downloads/Sequences/Alignments/"$s"_16KBph2.diploid.numts.fasta
done
done

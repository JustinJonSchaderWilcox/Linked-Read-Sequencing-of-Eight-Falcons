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
FalcoMito="/scratch/jw5478/FalconProject/Projects/RefGenome8/Mitochondria/Mitogenome"
HMMER="/home/jw5478/Executables/Programs/HMMER/bin"
Mafft="/home/jw5478/Executables/Programs/Mafft"
Emboss="/home/jw5478/Executables/Programs/EMBOSS-6.6.0/bin"

#Directories
rm -r $FalcoMito/Downloads/Sequences/Alignments/Consensus
mkdir $FalcoMito/Downloads/Sequences/Alignments/Consensus

for t in Peregrine Hierofalco Kestrel
do
#Date Variable
DATE=`ls $FalcoMito/Downloads/Sequences/"$t"MitoGenome_*.fasta | perl -pe "s#.*"$t"MitoGenome_(.*)\.fasta#\1#"`
#Index Downloaded Genomes
samtools faidx $FalcoMito/Downloads/Sequences/"$t"MitoGenome_"$DATE".fasta

#Alignment Directory
#Align for Consensus
$Mafft/mafft-linsi --nuc --adjustdirectionaccurately --reorder --thread 12 --maxiterate 1000 $FalcoMito/Downloads/Sequences/"$t"MitoGenome_"$DATE".fasta > $FalcoMito/Downloads/Sequences/Alignments/"$t"MitoGenome_"$DATE".linsi.fasta
awk 'NR==1 && />/{print $0; next}/>/{print "\n"$0}!/>/{printf $0}' $FalcoMito/Downloads/Sequences/"$t"MitoGenome_"$DATE".fasta | awk '/>/{print $0}!/>/{print $0$0}' > "$t"Mito/Downloads/Sequences/"$t"MitoGenome_"$DATE".concat.fasta

#HMMER Time
$HMMER/hmmbuild $FalcoMito/Downloads/Sequences/Alignments/"$t"MitoGenome_"$DATE".linsi.hmm  $FalcoMito/Downloads/Sequences/Alignments/"$t"MitoGenome_"$DATE".linsi.fasta
$HMMER/nhmmer --cpu 28 -A $FalcoMito/Downloads/Sequences/Alignments/"$t"MitoGenome_"$DATE".con_extract.sto -E 0.000001 --tblout $FalcoMito/Downloads/Sequences/Alignments/Records/"$t"MitoGenome_"$DATE".orient.txt --incE 0.000001   $FalcoMito/Downloads/Sequences/Alignments/"$t"MitoGenome_"$DATE".linsi.hmm $FalcoMito/Downloads/Sequences/"$t"MitoGenome_"$DATE".fasta

#Orient Sequences for Alignment
ACCESSIONS=`awk '!/^#/{print $1}' $FalcoMito/Downloads/Sequences/Alignments/Records/"$t"MitoGenome_"$DATE".orient.txt | sort | uniq`
rm $FalcoMito/Downloads/Sequences/"$t"MitoGenome_"$DATE".orient.fasta 
for a in $ACCESSIONS
do
DESCRIPTION=`awk -v accession=$a 'BEGIN{OFS="\t"}$1==accession{$1=$2=$3=$4=$5=$6=$7=$8=$9=$10=$11=$12=$13=$14=$15="";print $0; exit}' $FalcoMito/Downloads/Sequences/Alignments/Records/"$t"MitoGenome_"$DATE".orient.txt | perl -pe 's#\s+# #g'`
ORIENT=`awk -v accession=$a 'BEGIN{OFS="\t"}$1==accession{start=$5-1; exit}END{print start}' $FalcoMito/Downloads/Sequences/Alignments/Records/"$t"MitoGenome_"$DATE".orient.txt`
START=`expr $ORIENT + 1`
samtools faidx $FalcoMito/Downloads/Sequences/"$t"MitoGenome_"$DATE".fasta $a:$START > $FalcoMito/Downloads/Sequences/"$t"MitoGenome_"$DATE".orient_temp.fasta
samtools faidx $FalcoMito/Downloads/Sequences/"$t"MitoGenome_"$DATE".fasta $a:1-$ORIENT >> $FalcoMito/Downloads/Sequences/"$t"MitoGenome_"$DATE".orient_temp.fasta
awk -v accession=$a -v description="$DESCRIPTION" 'NR==1{print ">"accession" "description; next}!/>/{printf $0}END{printf "\n"}' $FalcoMito/Downloads/Sequences/"$t"MitoGenome_"$DATE".orient_temp.fasta >> $FalcoMito/Downloads/Sequences/"$t"MitoGenome_"$DATE".orient.fasta
done
rm $FalcoMito/Downloads/Sequences/"$t"MitoGenome_"$DATE".orient_temp.fasta
#Align Orient Circular Sequences and Produce Consensus Sequence
$Mafft/mafft-linsi --nuc --adjustdirectionaccurately --reorder --thread 12 --maxiterate 1000 $FalcoMito/Downloads/Sequences/"$t"MitoGenome_"$DATE".orient.fasta > $FalcoMito/Downloads/Sequences/Alignments/"$t"MitoGenome_"$DATE".orient.linsi.fasta
$Emboss/cons -sequence $FalcoMito/Downloads/Sequences/Alignments/"$t"MitoGenome_"$DATE".orient.linsi.fasta -outseq $FalcoMito/Downloads/Sequences/Alignments/Consensus/"$t"MitoGenome_"$DATE".cons.fasta
sed -i "s#>EMBOSS_001#>"$t"MitogenomeConsensus#" $FalcoMito/Downloads/Sequences/Alignments/Consensus/"$t"MitoGenome_"$DATE".cons.fasta
done

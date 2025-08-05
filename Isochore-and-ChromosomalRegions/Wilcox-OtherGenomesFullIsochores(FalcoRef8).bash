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

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Support
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Beds
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Seqs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Variants
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Variants/SNVs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Summary
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Summary/Composition

perl -pe 's#>.*chromosome ([0-9A-Z]+).*#>Chromosome_$1#' /scratch/jw5478/FalconProject/OtherGenomes/Catharus_ustulatus/GCF_009819885.1_bCatUst1.pri_genomic.fna | perl -pe 's#>.*(scaffold).*#>$1#' | perl -pe 's#>.*(mitochondion).*#>$1#' | awk 'BEGIN{OFS="."; count=0}$0==chrom{count=count+1; print chrom OFS count; next}/>Chromosome/{seq="INCLUDE"; chrom=$0; count=1; print $0 OFS count; next}/>scaffold/{seq="EXCLUDE"; next}/>mitochondrion/{seq="EXCLUDE"; next}seq=="INCLUDE"{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Thrush.chrom.fna &
perl -pe 's#>.*chromosome ([0-9A-Z]+).*#>Chromosome_$1#' /scratch/jw5478/FalconProject/OtherGenomes/Strigops_habroptila/GCF_004027225.2_bStrHab1.2.pri_genomic.fna | perl -pe 's#>.*(scaffold).*#>$1#' | perl -pe 's#>.*(mitochondion).*#>$1#' | awk 'BEGIN{OFS="."; count=0}$0==chrom{count=count+1; print chrom OFS count; next}/>Chromosome/{seq="INCLUDE"; chrom=$0; count=1; print $0 OFS count; next}/>scaffold/{seq="EXCLUDE"; next}/>mitochondrion/{seq="EXCLUDE"; next}seq=="INCLUDE"{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Kakapo.chrom.fna &
perl -pe 's#>.*chromosome ([0-9A-Z]+).*#>Chromosome_$1#' /scratch/jw5478/FalconProject/OtherGenomes/Cariama_cristata/GCA_009819825.1_bCarCri1.pri_genomic.fna | perl -pe 's#>.*(scaffold).*#>$1#' | perl -pe 's#>.*(mitochondion).*#>$1#' | awk 'BEGIN{OFS="."; count=0}$0==chrom{count=count+1; print chrom OFS count; next}/>Chromosome/{seq="INCLUDE"; chrom=$0; count=1; print $0 OFS count; next}/>scaffold/{seq="EXCLUDE"; next}/>mitochondrion/{seq="EXCLUDE"; next}seq=="INCLUDE"{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Seriema.chrom.fna &
perl -pe 's#>.*chromosome ([0-9A-Z]+).*#>Chromosome_$1#' /scratch/jw5478/FalconProject/OtherGenomes/Gallus_gallus/GCF_000002315.6_GRCg6a_genomic.fna | perl -pe 's#>.*(scaffold).*#>$1#' | perl -pe 's#>.*(mitochondion).*#>$1#' | awk 'BEGIN{OFS="."; count=0}$0==chrom{count=count+1; print chrom OFS count; next}/>Chromosome/{seq="INCLUDE"; chrom=$0; count=1; print $0 OFS count; next}/>scaffold/{seq="EXCLUDE"; next}/>mitochondrion/{seq="EXCLUDE"; next}seq=="INCLUDE"{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Gallus_gallus.chrom.fna &

perl -pe 's#>.*chromosome ([0-9A-Z]+).*#>Chromosome_$1#' /scratch/jw5478/FalconProject/OtherGenomes/Homo_sapiens/GCF_000001405.39_GRCh38.p13_genomic.fna | perl -pe 's#>.*(scaffold).*#>$1#' | perl -pe 's#>.*(mitochondion).*#>$1#' | awk 'BEGIN{OFS="."; count=0}$0==chrom{count=count+1; print chrom OFS count; next}/>Chromosome/{seq="INCLUDE"; chrom=$0; count=1; print $0 OFS count; next}/>scaffold/{seq="EXCLUDE"; next}/>mitochondrion/{seq="EXCLUDE"; next}seq=="INCLUDE"{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Human.chrom.fna &
perl -pe 's#>.*chromosome ([0-9A-Z]+).*#>Chromosome_$1#' /scratch/jw5478/FalconProject/OtherGenomes/Lacerta_agilus/GCF_009819535.1_rLacAgi1.pri_genomic.fna | perl -pe 's#>.*(scaffold).*#>$1#' | perl -pe 's#>.*(mitochondion).*#>$1#' | awk 'BEGIN{OFS="."; count=0}$0==chrom{count=count+1; print chrom OFS count; next}/>Chromosome/{seq="INCLUDE"; chrom=$0; count=1; print $0 OFS count; next}/>scaffold/{seq="EXCLUDE"; next}/>mitochondrion/{seq="EXCLUDE"; next}seq=="INCLUDE"{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Lacerta.chrom.fna &
perl -pe 's#>.*chromosome ([0-9A-Z]+).*#>Chromosome_$1#' /scratch/jw5478/FalconProject/OtherGenomes/Chelonia_mydas/GCF_015237465.1_rCheMyd1.pri_genomic.fna | perl -pe 's#>.*(scaffold).*#>$1#' | perl -pe 's#>.*(mitochondion).*#>$1#' | awk 'BEGIN{OFS="."; count=0}$0==chrom{count=count+1; print chrom OFS count; next}/>Chromosome/{seq="INCLUDE"; chrom=$0; count=1; print $0 OFS count; next}/>scaffold/{seq="EXCLUDE"; next}/>mitochondrion/{seq="EXCLUDE"; next}seq=="INCLUDE"{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Turtle.chrom.fna &
perl -pe 's#(>.*\.1).*#$1#' /scratch/jw5478/FalconProject/OtherGenomes/Crocodylas_porosis/GCF_001723895.1_CroPor_comp1_genomic.fna > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Crocodile.chrom.fna &
wait

samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Thrush.chrom.fna &
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Kakapo.chrom.fna &
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Seriema.chrom.fna &
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Gallus_gallus.chrom.fna &

samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Human.chrom.fna &
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Lacerta.chrom.fna &
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Turtle.chrom.fna &
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/Crocodile.chrom.fna &
wait

for o in Thrush Kakapo Seriema Gallus_gallus Human Lacerta Turtle Crocodile
do
awk 'BEGIN{OFS=FS="\t"}$2>=100000{print $1, $2}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/"$o".chrom.fna.fai | sort -k1,1V -k2,2n > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Support/"$o".chrom.genome
bedtools makewindows -g /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Support/"$o".chrom.genome -w 100000 | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2+1 OFS $3}' | awk '{Size=$3-$2}Size==99999{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Beds/"$o".chrom.win100KB.bed &
done
wait


#Function will extract sequences from these; paths need to be changed to extract from files for this analysis
function extract_regions {
WINDOW_POSITIONS=`cat /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Beds/"$o".chrom.win100KB.bed | perl -pe 's#\b\s+\b#__#g'`
rm /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Seqs/"$o".chrom.win100KB.fasta
for p in $WINDOW_POSITIONS
do
REF_SCAFFOLD=`echo $p | perl -pe 's#__#\t#g' | awk '{print $1}'`
REF_POS1=`echo $p | perl -pe 's#__#\t#g' | awk '{print $2}'`
REF_POS2=`echo $p | perl -pe 's#__#\t#g' | awk '{print $3}'`
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Genomes/"$o".chrom.fna $REF_SCAFFOLD:$REF_POS1-$REF_POS2 >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Seqs/"$o".chrom.win100KB.fasta
done
}

for o in Thrush Kakapo Seriema Gallus_gallus Human Lacerta Turtle Crocodile
do
extract_regions &
done
wait

#Pull out base composition of each Isochore
Seqtk="/home/jw5478/Executables/Programs/seqtk/seqtk"
for o in Thrush Kakapo Seriema Gallus_gallus Human Lacerta Turtle Crocodile
do
echo 'Win\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Summary/Composition/"$o".chrom.win100KB.comp
$Seqtk comp /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Seqs/"$o".chrom.win100KB.fasta >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Summary/Composition/"$o".chrom.win100KB.comp &
done
wait

function process_window { 
REF_SCAFFOLD=`echo $p | perl -pe 's#__#\t#g' | awk '{print $1}'`
REF_POS1=`echo $p | perl -pe 's#__#\t#g' | awk '{print $2}'`
REF_POS2=`echo $p | perl -pe 's#__#\t#g' | awk '{print $3}'`
WINDOW=`echo $p | perl -pe 's#__#:#' | perl -pe 's#__#-#'`

#Base Comp
A=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $3}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Summary/Composition/"$o".chrom.win100KB.comp`
C=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $4}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Summary/Composition/"$o".chrom.win100KB.comp`
G=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $5}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Summary/Composition/"$o".chrom.win100KB.comp`
T=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $6}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Summary/Composition/"$o".chrom.win100KB.comp`
CpG=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $9}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Summary/Composition/"$o".chrom.win100KB.comp`
tGC=`expr $G + $C`
tAT=`expr $A + $T`
tBASES=`expr $tAT + $tGC`
pGC=`echo -e "$tGC\t$tBASES" | awk '{print $1/$2}'`
echo -e "$o\t$p\t$REF_SCAFFOLD\t$COUNT\t$A\t$C\t$G\t$T\t$CpG\t$tGC\t$tBASES\t$pGC" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Summary/"$o"_WindowsSummary.txt
}

for o in Thrush Kakapo Seriema Gallus_gallus Human Lacerta Turtle Crocodile
do
echo -e "Genome\tWindow\tREF_SCAFFOLD\tWinNum\tA\tC\tG\tT\tCpG\ttGC\ttBASES\tpGC" > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Summary/"$o"_WindowsSummary.txt
WINDOW_POSITIONS=`cat /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/OtherGenomes_Full/Windows/Beds/"$o".chrom.win100KB.bed | perl -pe 's#\b\s+\b#__#g'`
COUNT=0
for p in $WINDOW_POSITIONS
do
process_window &
COUNT=`expr $COUNT + 1`
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 28 ]
do
sleep 1s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
done
wait

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Support
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Beds
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Seqs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Variants
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Variants/SNVs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Summary
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Summary/Composition

perl -pe 's#>.*chromosome ([0-9A-Z]+).*#>Chromosome_$1#' /scratch/jw5478/FalconProject/OtherGenomes/Catharus_ustulatus/GCF_009819885.1_bCatUst1.pri_genomic.fna | perl -pe 's#>.*(scaffold).*#>$1#' | perl -pe 's#>.*(mitochondion).*#>$1#' | awk 'BEGIN{OFS="."; count=0}$0==chrom{count=count+1; print chrom OFS count; next}/>Chromosome/{seq="INCLUDE"; chrom=$0; count=1; print $0 OFS count; next}/>scaffold/{seq="EXCLUDE"; next}/>mitochondrion/{seq="EXCLUDE"; next}seq=="INCLUDE"{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Thrush.chrom.fna &
perl -pe 's#>.*chromosome ([0-9A-Z]+).*#>Chromosome_$1#' /scratch/jw5478/FalconProject/OtherGenomes/Strigops_habroptila/GCF_004027225.2_bStrHab1.2.pri_genomic.fna | perl -pe 's#>.*(scaffold).*#>$1#' | perl -pe 's#>.*(mitochondion).*#>$1#' | awk 'BEGIN{OFS="."; count=0}$0==chrom{count=count+1; print chrom OFS count; next}/>Chromosome/{seq="INCLUDE"; chrom=$0; count=1; print $0 OFS count; next}/>scaffold/{seq="EXCLUDE"; next}/>mitochondrion/{seq="EXCLUDE"; next}seq=="INCLUDE"{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Kakapo.chrom.fna &
perl -pe 's#>.*chromosome ([0-9A-Z]+).*#>Chromosome_$1#' /scratch/jw5478/FalconProject/OtherGenomes/Cariama_cristata/GCA_009819825.1_bCarCri1.pri_genomic.fna | perl -pe 's#>.*(scaffold).*#>$1#' | perl -pe 's#>.*(mitochondion).*#>$1#' | awk 'BEGIN{OFS="."; count=0}$0==chrom{count=count+1; print chrom OFS count; next}/>Chromosome/{seq="INCLUDE"; chrom=$0; count=1; print $0 OFS count; next}/>scaffold/{seq="EXCLUDE"; next}/>mitochondrion/{seq="EXCLUDE"; next}seq=="INCLUDE"{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Seriema.chrom.fna &
perl -pe 's#>.*chromosome ([0-9A-Z]+).*#>Chromosome_$1#' /scratch/jw5478/FalconProject/OtherGenomes/Gallus_gallus/GCF_000002315.6_GRCg6a_genomic.fna | perl -pe 's#>.*(scaffold).*#>$1#' | perl -pe 's#>.*(mitochondion).*#>$1#' | awk 'BEGIN{OFS="."; count=0}$0==chrom{count=count+1; print chrom OFS count; next}/>Chromosome/{seq="INCLUDE"; chrom=$0; count=1; print $0 OFS count; next}/>scaffold/{seq="EXCLUDE"; next}/>mitochondrion/{seq="EXCLUDE"; next}seq=="INCLUDE"{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Gallus_gallus.chrom.fna &

perl -pe 's#>.*chromosome ([0-9A-Z]+).*#>Chromosome_$1#' /scratch/jw5478/FalconProject/OtherGenomes/Homo_sapiens/GCF_000001405.39_GRCh38.p13_genomic.fna | perl -pe 's#>.*(scaffold).*#>$1#' | perl -pe 's#>.*(mitochondion).*#>$1#' | awk 'BEGIN{OFS="."; count=0}$0==chrom{count=count+1; print chrom OFS count; next}/>Chromosome/{seq="INCLUDE"; chrom=$0; count=1; print $0 OFS count; next}/>scaffold/{seq="EXCLUDE"; next}/>mitochondrion/{seq="EXCLUDE"; next}seq=="INCLUDE"{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Human.chrom.fna &
perl -pe 's#>.*chromosome ([0-9A-Z]+).*#>Chromosome_$1#' /scratch/jw5478/FalconProject/OtherGenomes/Lacerta_agilus/GCF_009819535.1_rLacAgi1.pri_genomic.fna | perl -pe 's#>.*(scaffold).*#>$1#' | perl -pe 's#>.*(mitochondion).*#>$1#' | awk 'BEGIN{OFS="."; count=0}$0==chrom{count=count+1; print chrom OFS count; next}/>Chromosome/{seq="INCLUDE"; chrom=$0; count=1; print $0 OFS count; next}/>scaffold/{seq="EXCLUDE"; next}/>mitochondrion/{seq="EXCLUDE"; next}seq=="INCLUDE"{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Lacerta.chrom.fna &
perl -pe 's#>.*chromosome ([0-9A-Z]+).*#>Chromosome_$1#' /scratch/jw5478/FalconProject/OtherGenomes/Chelonia_mydas/GCF_015237465.1_rCheMyd1.pri_genomic.fna | perl -pe 's#>.*(scaffold).*#>$1#' | perl -pe 's#>.*(mitochondion).*#>$1#' | awk 'BEGIN{OFS="."; count=0}$0==chrom{count=count+1; print chrom OFS count; next}/>Chromosome/{seq="INCLUDE"; chrom=$0; count=1; print $0 OFS count; next}/>scaffold/{seq="EXCLUDE"; next}/>mitochondrion/{seq="EXCLUDE"; next}seq=="INCLUDE"{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Turtle.chrom.fna &
perl -pe 's#(>.*\.1).*#$1#' /scratch/jw5478/FalconProject/OtherGenomes/Crocodylas_porosis/GCF_001723895.1_CroPor_comp1_genomic.fna > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Crocodile.chrom.fna &
wait

samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Thrush.chrom.fna &
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Kakapo.chrom.fna &
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Seriema.chrom.fna &
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Gallus_gallus.chrom.fna &

samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Human.chrom.fna &
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Lacerta.chrom.fna &
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Turtle.chrom.fna &
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/Crocodile.chrom.fna &
wait

for o in Thrush Kakapo Seriema Gallus_gallus Human Lacerta Turtle Crocodile
do
awk 'BEGIN{OFS=FS="\t"}$2>=1000000{print $1, $2}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/"$o".chrom.fna.fai | sort -k1,1V -k2,2n > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Support/"$o".chrom.genome
bedtools makewindows -g /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Support/"$o".chrom.genome -w 1000000 | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2+1 OFS $3}' | awk '{Size=$3-$2}Size==999999{print $0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Beds/"$o".chrom.win100KB.bed &
done
wait


#Function will extract sequences from these; paths need to be changed to extract from files for this analysis
function extract_regions {
WINDOW_POSITIONS=`cat /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Beds/"$o".chrom.win100KB.bed | perl -pe 's#\b\s+\b#__#g'`
rm /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Seqs/"$o".chrom.win100KB.fasta
for p in $WINDOW_POSITIONS
do
REF_SCAFFOLD=`echo $p | perl -pe 's#__#\t#g' | awk '{print $1}'`
REF_POS1=`echo $p | perl -pe 's#__#\t#g' | awk '{print $2}'`
REF_POS2=`echo $p | perl -pe 's#__#\t#g' | awk '{print $3}'`
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Genomes/"$o".chrom.fna $REF_SCAFFOLD:$REF_POS1-$REF_POS2 >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Seqs/"$o".chrom.win100KB.fasta
done
}

for o in Thrush Kakapo Seriema Gallus_gallus Human Lacerta Turtle Crocodile
do
extract_regions &
done
wait

#Pull out base composition of each Isochore
Seqtk="/home/jw5478/Executables/Programs/seqtk/seqtk"
for o in Thrush Kakapo Seriema Gallus_gallus Human Lacerta Turtle Crocodile
do
echo 'Win\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Summary/Composition/"$o".chrom.win100KB.comp
$Seqtk comp /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Seqs/"$o".chrom.win100KB.fasta >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Summary/Composition/"$o".chrom.win100KB.comp &
done
wait

function process_window { 
REF_SCAFFOLD=`echo $p | perl -pe 's#__#\t#g' | awk '{print $1}'`
REF_POS1=`echo $p | perl -pe 's#__#\t#g' | awk '{print $2}'`
REF_POS2=`echo $p | perl -pe 's#__#\t#g' | awk '{print $3}'`
WINDOW=`echo $p | perl -pe 's#__#:#' | perl -pe 's#__#-#'`

#Base Comp
A=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $3}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Summary/Composition/"$o".chrom.win100KB.comp`
C=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $4}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Summary/Composition/"$o".chrom.win100KB.comp`
G=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $5}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Summary/Composition/"$o".chrom.win100KB.comp`
T=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $6}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Summary/Composition/"$o".chrom.win100KB.comp`
CpG=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $9}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Summary/Composition/"$o".chrom.win100KB.comp`
tGC=`expr $G + $C`
tAT=`expr $A + $T`
tBASES=`expr $tAT + $tGC`
pGC=`echo -e "$tGC\t$tBASES" | awk '{print $1/$2}'`
echo -e "$o\t$p\t$REF_SCAFFOLD\t$COUNT\t$A\t$C\t$G\t$T\t$CpG\t$tGC\t$tBASES\t$pGC" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Summary/"$o"_WindowsSummary.txt
}

for o in Thrush Kakapo Seriema Gallus_gallus Human Lacerta Turtle Crocodile
do
echo -e "Genome\tWindow\tREF_SCAFFOLD\tWinNum\tA\tC\tG\tT\tCpG\ttGC\ttBASES\tpGC" > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Summary/"$o"_WindowsSummary.txt
WINDOW_POSITIONS=`cat /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Megabases/OtherGenomes_Full/Windows/Beds/"$o".chrom.win100KB.bed | perl -pe 's#\b\s+\b#__#g'`
COUNT=0
for p in $WINDOW_POSITIONS
do
process_window &
COUNT=`expr $COUNT + 1`
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 28 ]
do
sleep 1s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
done
wait

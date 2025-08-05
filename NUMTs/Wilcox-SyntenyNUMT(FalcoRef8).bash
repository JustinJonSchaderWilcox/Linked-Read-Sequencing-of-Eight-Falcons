#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --mem=4GB
#SBATCH -p serial
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

module purge
module load perl
module load gencore/1
module load gencore_variant_detection/1.0

#Variables
FalconidaeMito="/scratch/jw5478/FalconProject/Projects/RefGenome8/Mitochondria/NUMT"

#Directory Structure
rm -r $FalconidaeMito/NUMT_Annotation
mkdir $FalconidaeMito/NUMT_Annotation

echo -e "Kestrel_ID\tGenome\tNUMT_SC\tScaffold_SC\tPos1_SC\tPos2_SC\tScaffold_Kest\tPos1_Kest\tPos2_Kest" > $FalconidaeMito/NUMT_Annotation/NUMT_Synteny.txt
for s in A B C D E F FpBSx77 Kestrel
do
NUMTs=`grep ">" $FalconidaeMito/Downloads/Sequences/Alignments/"$s"_16KBph2.diploid.numts.fasta | perl -pe 's#^>([0-9]+_[1-2]/[0-9]+-[0-9]+).*#$1#' | perl -pe 's#[/-]#\t#g' | awk 'BEGIN{OFS=FS="\t"}$3>$2{print $1 OFS $2 OFS $3}$3<$2{print $1 OFS $3 OFS $2}' | awk 'BEGIN{OFS=FS="\t"}$3>$2{size=$3-$2}$2>$3{size=$2-$3}size<14000{print $0}' | perl -pe 's#\t#__#g'`
for n in $NUMTs
do
SCAFF=`echo $n | perl -pe 's#__#\t#g' | awk '{print $1}'`
POS1=`echo $n | perl -pe 's#__#\t#g' | awk '{print $2}'`
POS2=`echo $n | perl -pe 's#__#\t#g' | awk '{print $3}'`
Kest_Coord=`awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 'BEGIN{OFS=FS="\t"}$15==scaff && $3>=pos1 && $13>0{kest_scaff=$14; kest_pos1=$1; kest_pos2$2}$15==scaff && $14==kest_scaff && $2>kest_pos2 && $4<pos2 && $13>0{kest_pos2=$2}$15==scaff && $14==kest_scaff && $4>=pos2 && kest_pos2<$2 && $13>0{kest_pos2=$2; exit}$15==scaff && $4>=pos1 && $13<0{kest_scaff=$14; kest_pos1=$1; kest_pos2=$2}$15==scaff && $14==kest_scaff && $2>kest_pos2 && $3<pos2 && $13<0{kest_pos2=$2}$15==scaff && $14==kest_scaff && $3>=pos2  && kest_pos2<$2 && $13<0{kest_pos2=$2; exit}END{print kest_scaff OFS kest_pos1 OFS kest_pos2}' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/MUMmerOutputs/"$s"-v-Kestrel100KBhap.coords`
Kest_Test=`echo -e "$Kest_Coord" | awk '{print $1}'` 
Kest_ID=`echo -e "$Kest_Coord" | perl -pe 's#\t#__#g'`
if [ -z $Kest_Test ]
then
Kest_Coord=`echo "Null\tNull\tNull"`
Kest_ID=Null
fi
echo -e "$Kest_ID\t$s\t$n\t$SCAFF\t$POS1\t$POS2\t$Kest_Coord" >> $FalconidaeMito/NUMT_Annotation/NUMT_Synteny.txt
done
done

for s in A B C D E F FpBSx77 Kestrel
do
awk -v genome=$s 'BEGIN{OFS=FS="\t"}$2==genome{print $7 OFS $8 OFS $9}'  $FalconidaeMito/NUMT_Annotation/NUMT_Synteny.txt | grep -v "Null" > $FalconidaeMito/NUMT_Annotation/"$s"_NUMT.bed
done
for s in A B C D E F FpBSx77 Kestrel
do
OTHER=`ls $FalconidaeMito/NUMT_Annotation/*_NUMT.bed | grep -v ""$s"_NUMT.bed"`
bedtools intersect -u -a $FalconidaeMito/NUMT_Annotation/"$s"_NUMT.bed -b $OTHER > $FalconidaeMito/NUMT_Annotation/"$s"_NUMT.shared.bed
bedtools subtract -A -a $FalconidaeMito/NUMT_Annotation/"$s"_NUMT.bed -b $OTHER > $FalconidaeMito/NUMT_Annotation/"$s"_NUMT.uniq.bed
done

echo -e "Genome\tTotalNUMTs\tNullNUMTs\tSharedNUMTs\tUniqNUMTs\thetNUMTs" > $FalconidaeMito/NUMT_Annotation/NUMT_SummaryReport.txt
for s in A B C D E F FpBSx77 Kestrel
do
TotalNUMTs=`awk -v genome=$s 'BEGIN{OFS=FS="\t"; count=0}$2==genome{count=count+1}END{print count}'  $FalconidaeMito/NUMT_Annotation/NUMT_Synteny.txt`
NullNUMTs=`awk -v genome=$s 'BEGIN{OFS=FS="\t"}$2==genome{print $7 OFS $8 OFS $9}'  $FalconidaeMito/NUMT_Annotation/NUMT_Synteny.txt | grep -c "Null"`
SharedNUMTs=`wc $FalconidaeMito/NUMT_Annotation/"$s"_NUMT.shared.bed | awk '{print $1}'`
UniqNUMTs=`wc $FalconidaeMito/NUMT_Annotation/"$s"_NUMT.uniq.bed | awk '{print $1}'`
homNUMTs=`awk -v genome=$s 'BEGIN{OFS=FS="\t"; count=0}$2==genome && $4~/1$/{allele=$1; next}$2==genome && $4~/2$/ && $1==allele{count=count+2}END{print count}' $FalconidaeMito/NUMT_Annotation/NUMT_Synteny.txt`
hetNUMTs=`expr $TotalNUMTs - $homNUMTs`
echo -e "$s\t$TotalNUMTs\t$NullNUMTs\t$SharedNUMTs\t$UniqNUMTs\t$hetNUMTs" >> $FalconidaeMito/NUMT_Annotation/NUMT_SummaryReport.txt
done
cat $FalconidaeMito/NUMT_Annotation/NUMT_SummaryReport.txt | perl -pe 's#^A\b#Lanner#' | perl -pe 's#^B\b#Barbary#' | perl -pe 's#^C\b#Peregrine#' | perl -pe 's#^D\b#Saker#' | perl -pe 's#^E\b#Gyr-1#' | perl -pe 's#^F\b#Gyr-2#' | perl -pe 's#^FpBSx77\b#Bl. Shaheen#' | perl -pe 's#NUMTs##g' | perl -pe 's#Null#Unaligned#' | perl -pe 's#Uniq#Unique#' > $FalconidaeMito/NUMT_Annotation/NUMT_UltSummaryReport.txt

echo -e "Genome\tTotalNUMTsBP\tNullNUMTsBP\tSharedNUMTsBP\tUniqNUMTsBP\thetNUMTsBP" > $FalconidaeMito/NUMT_Annotation/NUMT_SummaryReportBP.txt
for s in A B C D E F FpBSx77 Kestrel
do
TotalNUMTsBP=`awk -v genome=$s 'BEGIN{OFS=FS="\t";totalBP=0}$2==genome{size=$6-$5; totalBP=totalBP+size}END{print totalBP}'  $FalconidaeMito/NUMT_Annotation/NUMT_Synteny.txt`
NullNUMTsBP=`awk -v genome=$s 'BEGIN{OFS=FS="\t";totalBP=0}$2==genome && $1=="Null"{size=$6-$5; totalBP=totalBP+size}END{print totalBP}' $FalconidaeMito/NUMT_Annotation/NUMT_Synteny.txt`
homNUMTsBP=`awk -v genome=$s 'BEGIN{OFS=FS="\t"; homBP=0}$2==genome && $4~/1$/{allele=$1; size=$6-$5; next}$2==genome && $4~/2$/ && $1==allele{homBP=homBP+size; size=$6-$5; homBP=homBP+size}END{print homBP}' $FalconidaeMito/NUMT_Annotation/NUMT_Synteny.txt`
hetNUMTsBP=`expr $TotalNUMTsBP - $homNUMTsBP`

SharedNUMTsBP=0
for n in $(perl -pe 's#\t#__#g' $FalconidaeMito/NUMT_Annotation/"$s"_NUMT.shared.bed | sort -V | awk '$0!=prev{prev=$0; count=1; print $0"."count; next}$0==prev{count=count+1; print $0"."count}')
do
NUMT=`echo $n | perl -pe 's#\.#\t#' | awk '{print $1}'`
NUMBER=`echo $n | perl -pe 's#\.#\t#' | awk '{print $2}'`
tSharedNUMTsBP=`awk -v genome=$s -v numt=$NUMT -v number=$NUMBER 'BEGIN{OFS=FS="\t"; count=0}$2==genome && $1==numt{count=count+1}count==number{size=$6-$5; exit}END{print size}' $FalconidaeMito/NUMT_Annotation/NUMT_Synteny.txt`
SharedNUMTsBP=`expr $SharedNUMTsBP + $tSharedNUMTsBP`
done

UniqNUMTsBP=0
for n in $(perl -pe 's#\t#__#g' $FalconidaeMito/NUMT_Annotation/"$s"_NUMT.uniq.bed | sort -V | awk '$0!=prev{prev=$0; count=1; print $0"."count; next}$0==prev{count=count+1; print $0"."count}')
do
NUMT=`echo $n | perl -pe 's#\.#\t#' | awk '{print $1}'`
NUMBER=`echo $n | perl -pe 's#\.#\t#' | awk '{print $2}'`
tUniqNUMTsBP=`awk -v genome=$s -v numt=$NUMT -v number=$NUMBER 'BEGIN{OFS=FS="\t"; count=0}$2==genome && $1==numt{count=count+1}count==number{size=$6-$5; exit}END{print size}' $FalconidaeMito/NUMT_Annotation/NUMT_Synteny.txt`
UniqNUMTsBP=`expr $UniqNUMTsBP + $tUniqNUMTsBP`
done
echo -e "$s\t$TotalNUMTsBP\t$NullNUMTsBP\t$SharedNUMTsBP\t$UniqNUMTsBP\t$hetNUMTsBP" >> $FalconidaeMito/NUMT_Annotation/NUMT_SummaryReportBP.txt
done
cat $FalconidaeMito/NUMT_Annotation/NUMT_SummaryReportBP.txt | perl -pe 's#^A\b#Lanner#' | perl -pe 's#^B\b#Barbary#' | perl -pe 's#^C\b#Peregrine#' | perl -pe 's#^D\b#Saker#' | perl -pe 's#^E\b#Gyr-1#' | perl -pe 's#^F\b#Gyr-2#' | perl -pe 's#^FpBSx77\b#Bl. Shaheen#' | perl -pe 's#NUMTs##g' | perl -pe 's#Null#Unaligned#' | perl -pe 's#Uniq#Unique#' > $FalconidaeMito/NUMT_Annotation/NUMT_UltSummaryReportBP.txt

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
module load jdk

FalconidaeMito="/scratch/jw5478/FalconProject/Projects/RefGenome8/Mitochondria/NUMT"
HMMER="/home/jw5478/Executables/Programs/HMMER/bin"
Mafft="/home/jw5478/Executables/Programs/Mafft"
DATE=`ls $FalconidaeMito/Downloads/Sequences/AvesSampledFalconidaeMitoGenome_*.orient.fasta | perl -pe 's#.*FalconidaeMitoGenome_(.*)\.orient.fasta#\1#'`
FalcoMito="/scratch/jw5478/FalconProject/Projects/RefGenome8/Mitochondria/Mitogenome"

rm -r $FalconidaeMito/NUMT_Phylogeny/
mkdir $FalconidaeMito/NUMT_Phylogeny
mkdir $FalconidaeMito/NUMT_Phylogeny/Alignment
mkdir $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters
mkdir $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp
mkdir $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Aligned

$Mafft/mafft-linsi --nuc --adjustdirectionaccurately --reorder --thread 8 --maxiterate 1000 $FalconidaeMito/Downloads/Sequences/AvesSampledFalconidaeMitoGenome_"$DATE".orient.fasta > $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/AvesSampledFalconidaeMitoGenome_"$DATE".orient.linsi.fasta
$HMMER/hmmbuild $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/AvesSampledFalconidaeMitoGenome_"$DATE".orient.linsi.hmm  $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/AvesSampledFalconidaeMitoGenome_"$DATE".orient.linsi.fasta

for s in A B C D E F FpBSx77 Kestrel
do
perl -pe 's#(^>[0-9]+_[1-2]/[0-9]+-[0-9]+).*#$1#' $FalconidaeMito/Downloads/Sequences/Alignments/"$s"_16KBph2.diploid.numts.fasta | perl -pe "s#>#>$s:#" | perl -pe 's#[/-]#\t#g' | awk 'NR==1{print $0}NR!=1 && /^>/{print "\n"$0}!/>/{printf $0}' | awk 'BEGIN{OFS=FS="\t"}numt=="TRUE" && !/^>/{print $0; numt=""; next}$3>$2 && /^>/{size=$3-$2}$2>$3 && /^>/{size=$2-$3}size<14000 && /^>/{print $0; numt="TRUE"}' | perl -pe 's#\t#__#g' | perl -pe 's#_R_##' >> $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/AllFalco_16KBph2.diploid.numts.fasta
done
/home/jw5478/Executables/Programs/USEARCH/usearch11.0.667_i86linux32 -sortbylength $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/AllFalco_16KBph2.diploid.numts.fasta -fastaout $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/AllFalco_16KBph2.diploid.numts.sort_length.fasta
/home/jw5478/Executables/Programs/USEARCH/usearch11.0.667_i86linux32 -cluster_fast $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/AllFalco_16KBph2.diploid.numts.sort_length.fasta -threads 28 -id 0.8 -sort length -maxaccepts 1 -maxrejects 0 -strand both -sizeout -centroids $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/AllFalco_16KBph2.diploid.numts.centID80.fasta -uc $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/AllFalco_16KBph2.diploid.numts.centID80.uc
awk 'BEGIN{count=0}/^>/{print ">Cluster"count; count=count+1}!/>/{print $0}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/AllFalco_16KBph2.diploid.numts.centID80.fasta > $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/AllFalco_16KBph2.diploid.numts.centID80.clust.fasta
Clusters=`grep ">" $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/AllFalco_16KBph2.diploid.numts.centID80.clust.fasta | perl -pe 's#>##' | perl -pe 's#;.*##' | sort -V`
for f in $Clusters
do
Number=`echo $f | perl -pe 's#Cluster##'`
Centroid=`awk -v number="$Number" '{FS="\t"}$1=="S" && $10=="*" && $2==number{print $9; exit}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/AllFalco_16KBph2.diploid.numts.centID80.uc`
echo $Centroid > $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/TempCentHitsAllFalco_16KBph2.diploid.numts.txt
awk -v number="$Number" -v centroid="$Centroid" '{FS="\t"}$1=="H" && $10==centroid && $2==number{print $9}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/AllFalco_16KBph2.diploid.numts.centID80.uc >> $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/TempCentHitsAllFalco_16KBph2.diploid.numts.txt
/home/jw5478/Executables/Programs/USEARCH/usearch11.0.667_i86linux32 -fastx_getseqs $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/AllFalco_16KBph2.diploid.numts.fasta -labels $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/TempCentHitsAllFalco_16KBph2.diploid.numts.txt -fastaout $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/TempCentHitsAllFalco_16KBph2.diploid.numts.fasta
cat $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/TempCentHitsAllFalco_16KBph2.diploid.numts.fasta $FalconidaeMito/Downloads/Sequences/AvesSampledFalconidaeMitoGenome_"$DATE".orient.fasta $FalcoMito/Consensus/Mitogenomes/FalcoRef8_MitoGenome.cons.fasta > $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.fasta
#HMMER Time
$HMMER/nhmmer --cpu 12 -A $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.stk -E 0.000001 --tblout $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.orient.txt --incE 0.000001 $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/AvesSampledFalconidaeMitoGenome_"$DATE".orient.linsi.hmm $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.fasta
ACCESSIONS=`awk '!/^#/{print $1}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.orient.txt | sort | uniq`
rm $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.orient.fasta
for a in $ACCESSIONS
do
DESCRIPTION=`awk -v accession=$a 'BEGIN{OFS="\t"}$1==accession{$1=$2=$3=$4=$5=$6=$7=$8=$9=$10=$11=$12=$13=$14=$15="";print $0; exit}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.orient.txt | perl -pe 's#\s+# #g'`
ORIENT=`awk -v accession=$a 'BEGIN{OFS="\t"}$1==accession && $5>=$7{start=$5-$7; exit}$1==accession && $7>=$5{start=$7-$5; exit}END{print start}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.orient.txt`
START=`expr $ORIENT + 1`
samtools faidx $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.fasta $a:$START > $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.orient.temp.fasta
samtools faidx $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.fasta $a:1-$ORIENT >> $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.orient.temp.fasta
awk -v accession=$a -v description="$DESCRIPTION" 'NR==1{print ">"accession" "description; next}!/>/{printf $0}END{printf "\n"}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.orient.temp.fasta >> $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.orient.fasta
done
rm $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.orient.temp.fasta
$Mafft/mafft-linsi --nuc --adjustdirectionaccurately --reorder --thread 8 --maxiterate 1000 $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.orient.fasta | perl -pe 's#_R_##' > $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Aligned/"$f"AllFalco_16KBph2.diploid.numts.linsi.fasta
rm $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/TempCentHitsAllFalco_16KBph2.diploid.numts.fasta
rm $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp/"$f"AllFalco_16KBph2.diploid.numts.aves_ref8mitogenome.fasta
done
rm -r $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Temp

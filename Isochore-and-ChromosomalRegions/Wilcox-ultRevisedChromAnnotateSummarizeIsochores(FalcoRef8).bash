#!/bin/bash
#SBATCH --time=24:40:00
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

s=$1

function process_window { 
REF_SCAFFOLD=`echo $p | perl -pe 's#__#\t#g' | awk '{print $1}'`
REF_POS1=`echo $p | perl -pe 's#__#\t#g' | awk '{print $2}'`
REF_POS2=`echo $p | perl -pe 's#__#\t#g' | awk '{print $3}'`
WINDOW=`echo $p | perl -pe 's#__#:#' | perl -pe 's#__#-#'`

#Base Comp
A=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $3}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Summary/Composition/"$s"_16KBph2.haploid.win100KB.comp`
C=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $4}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Summary/Composition/"$s"_16KBph2.haploid.win100KB.comp`
G=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $5}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Summary/Composition/"$s"_16KBph2.haploid.win100KB.comp`
T=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $6}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Summary/Composition/"$s"_16KBph2.haploid.win100KB.comp`
CpG=`awk -v win=$WINDOW '{OFS=FS="\t"}NR==1{next}$1==win{print $9}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Summary/Composition/"$s"_16KBph2.haploid.win100KB.comp`
tGC=`expr $G + $C`
tAT=`expr $A + $T`
tBASES=`expr $tAT + $tGC`
pGC=`echo -e "$tGC\t$tBASES" | awk '{print $1/$2}'`

hetSNVs=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.norm.vcf.gz | grep -v "#" | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 'BEGIN{count=0}$1==scaff && $2>=pos1 && $2<=pos2 && length($4)==1 && length($5)==1 && $7=="PASS"{count=count+1}END{print count}'`
hetIndels=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.norm.vcf.gz | grep -v "#" | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 'BEGIN{count=0}$1==scaff && $2>=pos1 && $2<=pos2 && ((length($4)>1 && length($4)<=50 && length($5)==1) || (length($5)>1 && length($5)<=50 && length($4)==1)) && $7=="PASS"{count=count+1}END{print count}'`
hetIndelSVs=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.norm.vcf.gz /scratch/jw5478/FalconProject/Projects/RefGenome8/LongrangerOut/SelfAlign/16KBph2Haploid/"$s"_16KBph2-HaploidLR/outs/dels.vcf.gz | grep -v "#" | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 'BEGIN{count=0}$1==scaff && $2>=pos1 && $2<=pos2 && ((length($4)>1 && length($4)<=50 && length($5)==1) || (length($5)>1 && length($5)<=50 && length($4)==1)) && $7=="PASS"{count=count+1}END{print count}'`

#Unique Structural variation in each region using MUMmer alignments and size of insertions and deletions in each region using MUMmer alingment
uDelSV=`bedtools intersect -F 0.5 -a /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Beds/"$s"_16KBph2.haploid.win100KB.bed -b /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV."$s"coord.deletion50.uniq.bed | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$1==scaff && $2>=pos1 && $3<=pos2{print $0}' | wc | awk '{print $1}'`
uInsSV=`bedtools intersect -F 0.5 -a /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Beds/"$s"_16KBph2.haploid.win100KB.bed -b /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV."$s"coord.insertion50.uniq.bed | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$1==scaff && $2>=pos1 && $3<=pos2{print $0}' | wc | awk '{print $1}'`
uInvSV=`bedtools intersect -F 0.5 -a /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Beds/"$s"_16KBph2.haploid.win100KB.bed -b /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV."$s"coord.inversion.uniq.bed | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$1==scaff && $2>=pos1 && $3<=pos2{print $0}' | wc | awk '{print $1}'`
uSeg_DelSV=`bedtools intersect -F 0.5 -a /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Beds/"$s"_16KBph2.haploid.win100KB.bed -b /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV."$s"coord.seg_deletion.uniq.bed | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$1==scaff && $2>=pos1 && $3<=pos2{print $0}' | wc | awk '{print $1}'`
uSeg_InsSV=`bedtools intersect -F 0.5 -a /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Beds/"$s"_16KBph2.haploid.win100KB.bed -b /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV."$s"coord.seg_duplication.uniq.bed | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$1==scaff && $2>=pos1 && $3<=pos2{print $0}' | wc | awk '{print $1}'`

BPuInsSV=`bedtools intersect -F 0.5 -a /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Beds/"$s"_16KBph2.haploid.win100KB.bed -b /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV."$s"coord.insertion50.uniq.bed | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$1==scaff && $2>=pos1 && $3<=pos2{print $0}' | awk 'BEGIN{basepairs=0}{size=$3-$2; basepairs=basepairs+size}END{print basepairs}'`
BPuInvSV=`bedtools intersect -F 0.5 -a /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Beds/"$s"_16KBph2.haploid.win100KB.bed -b /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV."$s"coord.inversion.uniq.bed | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$1==scaff && $2>=pos1 && $3<=pos2{print $0}' | awk 'BEGIN{basepairs=0}{size=$3-$2; basepairs=basepairs+size}END{print basepairs}'`
BPuSeg_InsSV=`bedtools intersect -F 0.5 -a /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Beds/"$s"_16KBph2.haploid.win100KB.bed -b /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV."$s"coord.seg_duplication.uniq.bed | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$1==scaff && $2>=pos1 && $3<=pos2{print $0}' | awk 'BEGIN{basepairs=0}{size=$3-$2; basepairs=basepairs+size}END{print basepairs}'`

#Repeats in each region by type
rDNA=`grep -P "[_0-9+C]+" /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/"$s"/"$s"_16KBph2.diploid.fasta.out | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$5==scaff && $6>=pos1 && $7<=pos2{print $0}' | grep -c "DNA"`
rLTR=`grep -P "[_0-9+C]+" /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/"$s"/"$s"_16KBph2.diploid.fasta.out | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$5==scaff && $6>=pos1 && $7<=pos2{print $0}' | grep -c "LTR"`
rLINE=`grep -P "[_0-9+C]+" /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/"$s"/"$s"_16KBph2.diploid.fasta.out | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$5==scaff && $6>=pos1 && $7<=pos2{print $0}' | grep -c "LINE"`
rSINE=`grep -P "[_0-9+C]+" /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/"$s"/"$s"_16KBph2.diploid.fasta.out | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$5==scaff && $6>=pos1 && $7<=pos2{print $0}' | grep -c "SINE"`
rUnknown=`grep -P "[_0-9+C]+" /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/"$s"/"$s"_16KBph2.diploid.fasta.out | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$5==scaff && $6>=pos1 && $7<=pos2{print $0}' | grep -c "Unknown"`
BPrDNA=`grep -P "[_0-9+C]+" /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/"$s"/"$s"_16KBph2.diploid.fasta.out | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 'BEGIN{Bases=0; Count=0}$5==scaff && $6>=pos1 && $7<=pos2 && $11~/DNA/{Bases=$7-$6; Count=Count+Bases}END{print Count}'`
BPrLTR=`grep -P "[_0-9+C]+" /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/"$s"/"$s"_16KBph2.diploid.fasta.out | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 'BEGIN{Bases=0; Count=0}$5==scaff && $6>=pos1 && $7<=pos2 && $11~/LTR/{Bases=$7-$6; Count=Count+Bases}END{print Count}'`
BPrLINE=`grep -P "[_0-9+C]+" /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/"$s"/"$s"_16KBph2.diploid.fasta.out | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 'BEGIN{Bases=0; Count=0}$5==scaff && $6>=pos1 && $7<=pos2 && $11~/LINE/{Bases=$7-$6; Count=Count+Bases}END{print Count}'`
BPrSINE=`grep -P "[_0-9+C]+" /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/"$s"/"$s"_16KBph2.diploid.fasta.out | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 'BEGIN{Bases=0; Count=0}$5==scaff && $6>=pos1 && $7<=pos2 && $11~/SINE/{Bases=$7-$6; Count=Count+Bases}END{print Count}'`
BPrUnknown=`grep -P "[_0-9+C]+" /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/"$s"/"$s"_16KBph2.diploid.fasta.out | awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 'BEGIN{Bases=0; Count=0}$5==scaff && $6>=pos1 && $7<=pos2 && $11~/Unknown/{Bases=$7-$6; Count=Count+Bases}END{print Count}'`
#Full Length Gene Count
GENES=`awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$1==scaff && $2>=pos1 && $3<=pos2{print $0}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Transcriptome/Trintonate/A/BED/Transcripts/"$s"_tritonateTranscript.bed | wc | awk '{print $1}'`

#Chromosomal State
SCAFF_SIZE=`awk -v scaff=$REF_SCAFFOLD '$1==scaff{print $2}' /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta.fai`

if [[ $SCAFF_SIZE -ge 20000000 ]]
then
SIZE_CAT=`echo Large`
fi
if [[ $SCAFF_SIZE -lt 20000000 ]]
then
SIZE_CAT=`echo Small`
fi

THRUSH_CHROM_MATCH=`awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$15==scaff && $3>=pos1 && $4<=pos2 && $13>0{print $14}$15==scaff && $3<=pos2 && $4>=pos1 && $13<0{print $14}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-ThrushChrom.coords | sort | uniq | perl -pe 's#\n#__#' | perl -pe 's#__$##'`
KAKAPO_CHROM_MATCH=`awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$15==scaff && $3>=pos1 && $4<=pos2 && $13>0{print $14}$15==scaff && $3<=pos2 && $4>=pos1 && $13<0{print $14}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-KakapoChrom.coords | sort | uniq | perl -pe 's#\n#__#' | perl -pe 's#__$##'`
SERIEMA_CHROM_MATCH=`awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$15==scaff && $3>=pos1 && $4<=pos2 && $13>0{print $14}$15==scaff && $3<=pos2 && $4>=pos1 && $13<0{print $14}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-SeriemaChrom.coords | sort | uniq | perl -pe 's#\n#__#' | perl -pe 's#__$##'`

GYR_CHROM_MATCH=`awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$15==scaff && $3>=pos1 && $4<=pos2 && $13>0{print $14}$15==scaff && $3<=pos2 && $4>=pos1 && $13<0{print $14}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-GyrChrom.coords | sort | uniq | perl -pe 's#\n#__#' | perl -pe 's#__$##'`
PEREGRINE_CHROM_MATCH=`awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 '$15==scaff && $3>=pos1 && $4<=pos2 && $13>0{print $14}$15==scaff && $3<=pos2 && $4>=pos1 && $13<0{print $14}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/MUMmerOutputs/"$s"-v-pseudoPeregrineChrom.coords | sort | uniq | perl -pe 's#\n#__#' | perl -pe 's#__$##'`

THRUSH_CLASS_MATCH=""
for c in $(echo $THRUSH_CHROM_MATCH | perl -pe 's#__#\t#g')
do
THRUSH_CLASS_MATCH=`echo -e "$THRUSH_CLASS_MATCH""\n"$(awk -v scaff=$c '$1==scaff{print $3}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Support/Thrush.chrom.classification.txt) | perl -pe 's#^\n##' | sort | uniq`
done
THRUSH_CLASS_MATCH=`echo "$THRUSH_CLASS_MATCH" | perl -pe 's#\n#__#'  | perl -pe 's#__$##'`

KAKAPO_CLASS_MATCH=""
for c in $(echo $KAKAPO_CHROM_MATCH | perl -pe 's#__#\t#g')
do
KAKAPO_CLASS_MATCH=`echo -e "$KAKAPO_CLASS_MATCH""\n"$(awk -v scaff=$c '$1==scaff{print $3}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Support/Kakapo.chrom.classification.txt) | perl -pe 's#^\n##' | sort | uniq`
done
KAKAPO_CLASS_MATCH=`echo -e "$KAKAPO_CLASS_MATCH" | perl -pe 's#\n#__#'  | perl -pe 's#__$##'`


SERIEMA_CLASS_MATCH=""
for c in $(echo $SERIEMA_CHROM_MATCH | perl -pe 's#__#\t#g')
do
SERIEMA_CLASS_MATCH=`echo -e "$SERIEMA_CLASS_MATCH""\n"$(awk -v scaff=$c '$1==scaff{print $3}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Support/Seriema.chrom.classification.txt) | perl -pe 's#^\n##' | sort | uniq`
done
SERIEMA_CLASS_MATCH=`echo -e "$SERIEMA_CLASS_MATCH" | perl -pe 's#\n#__#'  | perl -pe 's#__$##'`


GYR_CLASS_MATCH=""
for c in $(echo $GYR_CHROM_MATCH | perl -pe 's#__#\t#g')
do
GYR_CLASS_MATCH=`echo -e "$GYR_CLASS_MATCH""\n"$(awk -v scaff=$c '$1==scaff{print $3}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Support/Gyr.chrom.classification.txt) | perl -pe 's#^\n##' | sort | uniq`
done
GYR_CLASS_MATCH=`echo -e "$GYR_CLASS_MATCH" | perl -pe 's#\n#__#'  | perl -pe 's#__$##'`


PEREGRINE_CLASS_MATCH=""
for c in $(echo $PEREGRINE_CHROM_MATCH | perl -pe 's#__#\t#g')
do
PEREGRINE_CLASS_MATCH=`echo -e "$PEREGRINE_CLASS_MATCH""\n"$(awk -v scaff=$c '$1==scaff{print $3}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Support/Peregrine.chrom.classification.txt) | perl -pe 's#^\n##' | sort | uniq`
done
PEREGRINE_CLASS_MATCH=`echo -e "$PEREGRINE_CLASS_MATCH" | perl -pe 's#\n#__#'  | perl -pe 's#__$##'`


if [[ -z $THRUSH_CHROM_MATCH ]]
then
THRUSH_CHROM_MATCH=`echo NULL`
THRUSH_CLASS_MATCH=`echo NULL`
fi

if [[ -z $KAKAPO_CHROM_MATCH ]]
then
KAKAPO_CHROM_MATCH=`echo NULL`
KAKAPO_CLASS_MATCH=`echo NULL`
fi

if [[ -z $SERIEMA_CHROM_MATCH ]]
then
SERIEMA_CHROM_MATCH=`echo NULL`
SERIEMA_CLASS_MATCH=`echo NULL`
fi

OVERALL_CLASS_MATCH=`echo -e "$THRUSH_CLASS_MATCH\n$KAKAPO_CLASS_MATCH\n$SERIEMA_CLASS_MATCH" | perl -pe 's#__#\n#g' | sort | uniq | perl -pe 's#\n#__#'  | perl -pe 's#__$##' | perl -pe 's#(_)*NULL(_)*##'`

CLASS_MATCH_NUMBER=`echo -e "$THRUSH_CLASS_MATCH\n$KAKAPO_CLASS_MATCH\n$SERIEMA_CLASS_MATCH"  | perl -pe 's#__#\n#g' | sort | uniq | perl -pe 's#(\n)*NULL(\n)*##' | wc | awk '{print $1}'`

if [[ $CLASS_MATCH_NUMBER -gt 1 ]]
then
CONSENSUS_CLASS_MATCH=`echo Ambiguous`
fi
if [[ $CLASS_MATCH_NUMBER -eq 1 ]]
then
CONSENSUS_CLASS_MATCH=`echo $OVERALL_CLASS_MATCH`
fi

if [[ -z $GYR_CHROM_MATCH ]]
then
GYR_CHROM_MATCH=`echo NULL`
GYR_CLASS_MATCH=`echo NULL`
fi

if [[ -z $PEREGRINE_CHROM_MATCH ]]
then
PEREGRINE_CHROM_MATCH=`echo NULL`
PEREGRINE_CLASS_MATCH=`echo NULL`
fi

FALCO_CLASS_MATCH=`echo -e "$GYR_CLASS_MATCH\n$PEREGRINE_CLASS_MATCH" | perl -pe 's#__#\n#g' | sort | uniq | perl -pe 's#\n#__#'  | perl -pe 's#__$##' | perl -pe 's#(_)*NULL(_)*##'`
FalcoCLASS_MATCH_NUMBER=`echo -e "$GYR_CLASS_MATCH\n$PEREGRINE_CLASS_MATCH"  | perl -pe 's#__#\n#g' | sort | uniq | perl -pe 's#(\n)*NULL(\n)*##' | wc | awk '{print $1}'`

if [[ $FalcoCLASS_MATCH_NUMBER -gt 1 ]]
then
FalcoCONSENSUS_CLASS_MATCH=`echo Ambiguous`
fi
if [[ $FalcoCLASS_MATCH_NUMBER -eq 1 ]]
then
FalcoCONSENSUS_CLASS_MATCH=`echo $FALCO_CLASS_MATCH`
fi

FalcoSex_Chrom_Match=`echo $GYR_CHROM_MATCH $PEREGRINE_CHROM_MATCH | perl -pe 's#\s+#__#g' | perl -pe 's#__#\n#g' | awk '/Z/{print $0}/W/{print$0}' | sort | uniq | perl -pe 's#\n#__#g' | perl -pe 's#__$##'`
Num_Chrom_Match=`echo $GYR_CHROM_MATCH $PEREGRINE_CHROM_MATCH | perl -pe 's#\s+#__#g' | perl -pe 's#__#\n#g' | sort | uniq | wc | awk '{print $1}'`
Num_SexChrom_Match=`echo $FalcoSex_Chrom_Match | perl -pe 's#__#\n#g' | sort | uniq | wc | awk '{print $1}'`
NumExcessAutosome=`expr $Num_Chrom_Match - $Num_SexChrom_Match`

if [[ $Num_SexChrom_Match -gt 1 ]]
then
ConFalcoSex_Chrom_Type=`echo Ambiguous`
fi
if [[ $Num_SexChrom_Match -eq 1 ]]
then
ConFalcoSex_Chrom_Type=`echo $FalcoSex_Chrom_Match`
fi

if [[ -z $ConFalcoSex_Chrom_Type ]]
then
ConFalcoSex_Chrom_Type=`echo Autosome`
fi

if [[ $NumExcessAutosome -gt 0 ]]
then
ConFalcoSex_Chrom_Type=`echo Ambiguous`
fi


#SNPs
KESTREL_SCAFFOLD=`awk -v scaff=$REF_SCAFFOLD -v pos1=$REF_POS1 -v pos2=$REF_POS2 'BEGIN{OFS=FS="\t"}$15==scaff && $3>=pos1 && $4<=pos2 && $13>0{offset=$3-pos1; kest_pos1=$1+offset; kest_pos2=$2+offset; print $14 OFS kest_pos1 OFS kest_pos2}$15==scaff && $3<=pos2 && $4>=pos1 && $13<0{offset=$3-pos1; kest_pos1=$1+offset; kest_pos2=$2+offset; print $14 OFS kest_pos1 OFS kest_pos2}' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/MUMmerOutputs/"$s"-v-Kestrel100KBhap.coords | perl -pe 's#\t#__#g' | perl -pe 's#\n#___#' | perl -pe 's#___$##'`
t_wCpG_AT=0
t_wCpG_AC=0
t_wCpG_AG=0
t_wCpG_TA=0
t_wCpG_TC=0
t_wCpG_TG=0
t_wCpG_CA=0
t_wCpG_CT=0
t_wCpG_CG=0
t_wCpG_GA=0
t_wCpG_GT=0
t_wCpG_GC=0
t_nCpG_AT=0
t_nCpG_AC=0
t_nCpG_AG=0
t_nCpG_TA=0
t_nCpG_TC=0
t_nCpG_TG=0
t_nCpG_CA=0
t_nCpG_CT=0
t_nCpG_CG=0
t_nCpG_GA=0
t_nCpG_GT=0
t_nCpG_GC=0

f_wCpG_AT=0
f_wCpG_AC=0
f_wCpG_AG=0
f_wCpG_TA=0
f_wCpG_TC=0
f_wCpG_TG=0
f_wCpG_CA=0
f_wCpG_CT=0
f_wCpG_CG=0
f_wCpG_GA=0
f_wCpG_GT=0
f_wCpG_GC=0
f_nCpG_AT=0
f_nCpG_AC=0
f_nCpG_AG=0
f_nCpG_TA=0
f_nCpG_TC=0
f_nCpG_TG=0
f_nCpG_CA=0
f_nCpG_CT=0
f_nCpG_CG=0
f_nCpG_GA=0
f_nCpG_GT=0
f_nCpG_GC=0

wCPG_CombinedAT=0
wCPG_CombinedGC=0
nCPG_CombinedAT=0
nCPG_CombinedGC=0

for c in $(echo $KESTREL_SCAFFOLD | perl -pe 's#___#\t#g')
do
SCAFF=`echo $c | perl -pe 's#__#\t#g' | awk '{print $1}'`
POS1=`echo $c | perl -pe 's#__#\t#g' | awk '{print $2}'`
POS2=`echo $c | perl -pe 's#__#\t#g' | awk '{print $3}'`

BPuDelSV=`awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 'BEGIN{OFS=FS="\t"}$1==scaff && $2>=pos1 && $3<=pos2{print $0}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV.kestCoord.deletion50.uniq.bed | awk 'BEGIN{basepairs=0}{size=$3-$2; basepairs=basepairs+size}END{print basepairs}'`
BPuSeg_DelSV=`awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 'BEGIN{OFS=FS="\t"}$1==scaff && $2>=pos1 && $3<=pos2{print $0}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV.kestCoord.seg_duplication.uniq.bed | awk 'BEGIN{basepairs=0}{size=$3-$2; basepairs=basepairs+size}END{print basepairs}'`

###Base Composition
#Finding Unique Variant File for Sample
VCF=`grep "for stripped" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/keyREADME.txt | grep "$s-v-Kestrel" | perl -pe 's#.*([0-9]{4}.vcf).*#$1#'`
#'A>T'
t_wCpG_AT=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="A" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'A>C'
t_wCpG_AC=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="A" && $5=="C"{print $0}' | wc | awk '{print $1}'`
#'A>G'
t_wCpG_AG=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="A" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'T>A'
t_wCpG_TA=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="T" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'T>C'
t_wCpG_TC=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="T" && $5=="C"{print $0}' | wc | awk '{print $1}'`
#'T>G'
t_wCpG_TG=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="T" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'C>A'
t_wCpG_CA=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="C" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'C>T'
t_wCpG_CT=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="C" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'C>G'
t_wCpG_CG=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="C" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'G>A'
t_wCpG_GA=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="G" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'G>T'
t_wCpG_GT=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="G" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'G>C'
t_wCpG_GC=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="G" && $5=="C"{print $0}' | wc | awk '{print $1}'`

#Indels
tNumDeletion=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$1==scaff && $2>=pos1 && $2<=pos2 && $8=="DELETION"{print $0}' | wc | awk '{print $1}'`
tNumInsertion=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$1==scaff && $2>=pos1 && $2<=pos2 && $8=="INSERTION"{print $0}' | wc | awk '{print $1}'`
tDeletionBP=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 'BEGIN{Count=0}$1==scaff && $2>=pos1 && $2<=pos2 && $8=="DELETION"{Count=Count+length($4)}END{print Count}'`
tInsertionBP=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 'BEGIN{Count=0}$1==scaff && $2>=pos1 && $2<=pos2 && $8=="INSERTION"{Count=Count+length($5)}END{print Count}'`

fNumDeletion=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$1==scaff && $2>=pos1 && $2<=pos2 && $8=="DELETION"{print $0}' | wc | awk '{print $1}'`
fNumInsertion=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$1==scaff && $2>=pos1 && $2<=pos2 && $8=="INSERTION"{print $0}' | wc | awk '{print $1}'`
fDeletionBP=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 'BEGIN{Count=0}$1==scaff && $2>=pos1 && $2<=pos2 && $8=="DELETION"{Count=Count+length($4)}END{print Count}'`
fInsertionBP=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 'BEGIN{Count=0}$1==scaff && $2>=pos1 && $2<=pos2 && $8=="INSERTION"{Count=Count+length($5)}END{print Count}'`

##noCpG
#Finding Unique Variant File for Sample
#'A>T'
t_nCpG_AT=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="A" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'A>C'
t_nCpG_AC=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="A" && $5=="C"{print $0}' | wc | awk '{print $1}'`
#'A>G'
t_nCpG_AG=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="A" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'T>A'
t_nCpG_TA=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="T" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'T>C'
t_nCpG_TC=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="T" && $5=="C"{print $0}' | wc | awk '{print $1}'`
#'T>G'
t_nCpG_TG=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="T" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'C>A'
t_nCpG_CA=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="C" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'C>T'
t_nCpG_CT=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="C" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'C>G'
t_nCpG_CG=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="C" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'G>A'
t_nCpG_GA=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="G" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'G>T'
t_nCpG_GT=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="G" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'G>C'
t_nCpG_GC=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="G" && $5=="C"{print $0}' | wc | awk '{print $1}'`

###Base Composition
#Finding Unique Variant File for Sample
VCF=`grep "for stripped" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/keyREADME.txt | grep "$s-v-Kestrel" | perl -pe 's#.*([0-9]{4}.vcf).*#$1#'`
#'A>T'
f_wCpG_AT=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="A" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'A>C'
f_wCpG_AC=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="A" && $5=="C"{print $0}' | wc | awk '{print $1}'`
#'A>G'
f_wCpG_AG=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="A" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'T>A'
f_wCpG_TA=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="T" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'T>C'
f_wCpG_TC=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="T" && $5=="C"{print $0}' | wc | awk '{print $1}'`
#'T>G'
f_wCpG_TG=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="T" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'C>A'
f_wCpG_CA=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="C" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'C>T'
f_wCpG_CT=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="C" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'C>G'
f_wCpG_CG=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="C" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'G>A'
f_wCpG_GA=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="G" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'G>T'
f_wCpG_GT=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="G" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'G>C'
f_wCpG_GC=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="G" && $5=="C"{print $0}' | wc | awk '{print $1}'`

#Indels
NumDeletion=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$1==scaff && $2>=pos1 && $2<=pos2 && $8=="DELETION"{print $0}' | wc | awk '{print $1}'`
NumInsertion=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$1==scaff && $2>=pos1 && $2<=pos2 && $8=="INSERTION"{print $0}' | wc | awk '{print $1}'`
DeletionBP=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 'BEGIN{Count=0}$1==scaff && $2>=pos1 && $2<=pos2 && $8=="DELETION"{Count=Count+length($4)}END{print Count}'`
InsertionBP=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/$VCF | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 'BEGIN{Count=0}$1==scaff && $2>=pos1 && $2<=pos2 && $8=="INSERTION"{Count=Count+length($5)}END{print Count}'`

##noCpG
#Finding Unique Variant File for Sample
#'A>T'
f_nCpG_AT=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="A" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'A>C'
f_nCpG_AC=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="A" && $5=="C"{print $0}' | wc | awk '{print $1}'`
#'A>G'
f_nCpG_AG=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="A" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'T>A'
f_nCpG_TA=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="T" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'T>C'
f_nCpG_TC=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="T" && $5=="C"{print $0}' | wc | awk '{print $1}'`
#'T>G'
f_nCpG_TG=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="T" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'C>A'
f_nCpG_CA=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="C" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'C>T'
f_nCpG_CT=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="C" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'C>G'
f_nCpG_CG=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="C" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'G>A'
f_nCpG_GA=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="G" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'G>T'
f_nCpG_GT=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="G" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'G>C'
f_nCpG_GC=`grep "^$SCAFF" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/RevisedUniques/noGC/SNVs/"$s"/0000.vcf | awk '$0==prev{print $0}{prev=$0}' | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 '$2<pos1{next}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && $4=="G" && $5=="C"{print $0}' | wc | awk '{print $1}'`

t_wCpG_AT=`expr $t_wCpG_AT + $t_wCpG_AT`
t_wCpG_AC=`expr $t_wCpG_AC + $t_wCpG_AC`
t_wCpG_AG=`expr $t_wCpG_AG + $t_wCpG_AG`
t_wCpG_TA=`expr $t_wCpG_TA + $t_wCpG_TA`
t_wCpG_TC=`expr $t_wCpG_TC + $t_wCpG_TC`
t_wCpG_TG=`expr $t_wCpG_TG + $t_wCpG_TG`
t_wCpG_CA=`expr $t_wCpG_CA + $t_wCpG_CA`
t_wCpG_CT=`expr $t_wCpG_CT + $t_wCpG_CT`
t_wCpG_CG=`expr $t_wCpG_CG + $t_wCpG_CG`
t_wCpG_GA=`expr $t_wCpG_GA + $t_wCpG_GA`
t_wCpG_GT=`expr $t_wCpG_GT + $t_wCpG_GT`
t_wCpG_GC=`expr $t_wCpG_GC + $t_wCpG_GC`
t_nCpG_AT=`expr $t_nCpG_AT + $t_nCpG_AT`
t_nCpG_AC=`expr $t_nCpG_AC + $t_nCpG_AC`
t_nCpG_AG=`expr $t_nCpG_AG + $t_nCpG_AG`
t_nCpG_TA=`expr $t_nCpG_TA + $t_nCpG_TA`
t_nCpG_TC=`expr $t_nCpG_TC + $t_nCpG_TC`
t_nCpG_TG=`expr $t_nCpG_TG + $t_nCpG_TG`
t_nCpG_CA=`expr $t_nCpG_CA + $t_nCpG_CA`
t_nCpG_CT=`expr $t_nCpG_CT + $t_nCpG_CT`
t_nCpG_CG=`expr $t_nCpG_CG + $t_nCpG_CG`
t_nCpG_GA=`expr $t_nCpG_GA + $t_nCpG_GA`
t_nCpG_GT=`expr $t_nCpG_GT + $t_nCpG_GT`
t_nCpG_GC=`expr $t_nCpG_GC + $t_nCpG_GC`

f_wCpG_AT=`expr $f_wCpG_AT + $f_wCpG_AT`
f_wCpG_AC=`expr $f_wCpG_AC + $f_wCpG_AC`
f_wCpG_AG=`expr $f_wCpG_AG + $f_wCpG_AG`
f_wCpG_TA=`expr $f_wCpG_TA + $f_wCpG_TA`
f_wCpG_TC=`expr $f_wCpG_TC + $f_wCpG_TC`
f_wCpG_TG=`expr $f_wCpG_TG + $f_wCpG_TG`
f_wCpG_CA=`expr $f_wCpG_CA + $f_wCpG_CA`
f_wCpG_CT=`expr $f_wCpG_CT + $f_wCpG_CT`
f_wCpG_CG=`expr $f_wCpG_CG + $f_wCpG_CG`
f_wCpG_GA=`expr $f_wCpG_GA + $f_wCpG_GA`
f_wCpG_GT=`expr $f_wCpG_GT + $f_wCpG_GT`
f_wCpG_GC=`expr $f_wCpG_GC + $f_wCpG_GC`
f_nCpG_AT=`expr $f_nCpG_AT + $f_nCpG_AT`
f_nCpG_AC=`expr $f_nCpG_AC + $f_nCpG_AC`
f_nCpG_AG=`expr $f_nCpG_AG + $f_nCpG_AG`
f_nCpG_TA=`expr $f_nCpG_TA + $f_nCpG_TA`
f_nCpG_TC=`expr $f_nCpG_TC + $f_nCpG_TC`
f_nCpG_TG=`expr $f_nCpG_TG + $f_nCpG_TG`
f_nCpG_CA=`expr $f_nCpG_CA + $f_nCpG_CA`
f_nCpG_CT=`expr $f_nCpG_CT + $f_nCpG_CT`
f_nCpG_CG=`expr $f_nCpG_CG + $f_nCpG_CG`
f_nCpG_GA=`expr $f_nCpG_GA + $f_nCpG_GA`
f_nCpG_GT=`expr $f_nCpG_GT + $f_nCpG_GT`
f_nCpG_GC=`expr $f_nCpG_GC + $f_nCpG_GC`

tNumDeletion=`expr $tNumDeletion + $tNumDeletion`
fNumDeletion=`expr $fNumDeletion + $fNumDeletion`
tNumInsertion=`expr $tNumInsertion - $tNumInsertion`
fNumInsertion=`expr $fNumInsertion - $fNumInsertion`
tDeletionBP=`expr $tDeletionBP - $tDeletionBP`
fDeletionBP=`expr $fDeletionBP - $fDeletionBP`
tInsertionBP=`expr $tInsertionBP - $tInsertionBP`
fInsertionBP=`expr $fInsertionBP - $fInsertionBP`

LINE2_wCpG=`head -$POS2 /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Temp/gVCF_SCAFF/"$SCAFF"_wGC.pseudo_vcf | awk -v pos2=$POS2 '$2>=pos2{print NR; exit}'`
LINE1_nCpG=`head -$POS2 /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Temp/gVCF_SCAFF/"$SCAFF"_noGC.pseudo_vcf | awk -v pos1=$POS1 '$2>=pos1{print NR; exit}'`
DIFF_TOTAL=`expr $LINE2_wCpG - $LINE1_nCpG`

t_wCPG_CombinedAT=`head -$LINE2_wCpG /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Temp/gVCF_SCAFF/"$SCAFF"_wGC.pseudo_vcf | tail -$DIFF_TOTAL | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 'BEGIN{Count=0}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && ($4=="A" || $4=="T"){Count=Count+1}END{print Count}'`
t_wCPG_CombinedGC=`head -$LINE2_wCpG /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Temp/gVCF_SCAFF/"$SCAFF"_wGC.pseudo_vcf | tail -$DIFF_TOTAL | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 'BEGIN{Count=0}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && ($4=="G" || $4=="C"){Count=Count+1}END{print Count}'`
wCPG_CombinedAT=`expr $wCPG_CombinedAT + $t_wCPG_CombinedAT`
wCPG_CombinedGC=`expr $wCPG_CombinedGC + $t_wCPG_CombinedGC`

t_nCpG_CombinedAT=`head -$LINE2_wCpG /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Temp/gVCF_SCAFF/"$SCAFF"_noGC.pseudo_vcf | tail -$DIFF_TOTAL | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 'BEGIN{Count=0}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && ($4=="A" || $4=="T"){Count=Count+1}END{print Count}'`
t_nCpG_CombinedGC=`head -$LINE2_wCpG /scratch/jw5478/FalconProject/Projects/RefGenome8/Chromosomal/Chrom_Alignments/Temp/gVCF_SCAFF/"$SCAFF"_noGC.pseudo_vcf | tail -$DIFF_TOTAL | awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 'BEGIN{Count=0}$2>pos2{exit}$1==scaff && $2>=pos1 && $2<=pos2 && ($4=="G" || $4=="C"){Count=Count+1}END{print Count}'`
nCpG_CombinedAT=`expr $nCpG_CombinedAT + $t_nCpG_CombinedAT`
nCpG_CombinedGC=`expr $nCpG_CombinedGC + $t_nCpG_CombinedGC`
ROH=`awk -v scaff=$SCAFF -v pos1=$POS1 -v pos2=$POS2 'BEGIN{count=0}$1==scaff && $2>=pos1 && $2<=pos2 && $3==1 && $4>95{count=count+1}END{print count}' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/ROH/"$s"_16KBph2-Haploid.pv.roh`
done

h_wCpG_AT=`expr $t_wCpG_AT - $f_wCpG_AT`
h_wCpG_AC=`expr $t_wCpG_AC - $f_wCpG_AC`
h_wCpG_AG=`expr $t_wCpG_AG - $f_wCpG_AG`
h_wCpG_TA=`expr $t_wCpG_TA - $f_wCpG_TA`
h_wCpG_TC=`expr $t_wCpG_TC - $f_wCpG_TC`
h_wCpG_TG=`expr $t_wCpG_TG - $f_wCpG_TG`
h_wCpG_CA=`expr $t_wCpG_CA - $f_wCpG_CA`
h_wCpG_CT=`expr $t_wCpG_CT - $f_wCpG_CT`
h_wCpG_CG=`expr $t_wCpG_CG - $f_wCpG_CG`
h_wCpG_GA=`expr $t_wCpG_GA - $f_wCpG_GA`
h_wCpG_GT=`expr $t_wCpG_GT - $f_wCpG_GT`
h_wCpG_GC=`expr $t_wCpG_GC - $f_wCpG_GC`

h_nCpG_AT=`expr $t_nCpG_AT - $f_nCpG_AT`
h_nCpG_AC=`expr $t_nCpG_AC - $f_nCpG_AC`
h_nCpG_AG=`expr $t_nCpG_AG - $f_nCpG_AG`
h_nCpG_TA=`expr $t_nCpG_TA - $f_nCpG_TA`
h_nCpG_TC=`expr $t_nCpG_TC - $f_nCpG_TC`
h_nCpG_TG=`expr $t_nCpG_TG - $f_nCpG_TG`
h_nCpG_CA=`expr $t_nCpG_CA - $f_nCpG_CA`
h_nCpG_CT=`expr $t_nCpG_CT - $f_nCpG_CT`
h_nCpG_CG=`expr $t_nCpG_CG - $f_nCpG_CG`
h_nCpG_GA=`expr $t_nCpG_GA - $f_nCpG_GA`
h_nCpG_GT=`expr $t_nCpG_GT - $f_nCpG_GT`
h_nCpG_GC=`expr $t_nCpG_GC - $f_nCpG_GC`

hNumDeletion=`expr $tNumDeletion - $fNumDeletion`
hNumInsertion=`expr $tNumInsertion - $fNumInsertion`
hDeletionBP=`expr $tDeletionBP - $fDeletionBP`
hInsertionBP=`expr $tInsertionBP - $fInsertionBP`

#Summarize
echo -e "$s\t$p\t$REF_SCAFFOLD\t$COUNT\t$A\t$C\t$G\t$T\t$CpG\t$tGC\t$tBASES\t$pGC\t$uDelSV\t$uInsSV\t$uInvSV\t$uSeg_DelSV\t$uSeg_InsSV\t$BPuDelSV\t$BPuInsSV\t$BPuInvSV\t$BPuSeg_DelSV\t$BPuSeg_InsSV\t$rDNA\t$rLTR\t$rLINE\t$rSINE\t$rUnknown\t$BPrDNA\t$BPrLTR\t$BPrLINE\t$BPrSINE\t$BPrUnknown\t$GENES\t$SCAFF_SIZE\t$SIZE_CAT\t$THRUSH_CHROM_MATCH\t$KAKAPO_CHROM_MATCH\t$SERIEMA_CHROM_MATCH\t$THRUSH_CLASS_MATCH\t$KAKAPO_CLASS_MATCH\t$SERIEMA_CLASS_MATCH\t$OVERALL_CLASS_MATCH\t$CLASS_MATCH_NUMBER\t$CONSENSUS_CLASS_MATCH\t$f_wCpG_AT\t$f_wCpG_AG\t$f_wCpG_AC\t$f_wCpG_TA\t$f_wCpG_TC\t$f_wCpG_TG\t$f_wCpG_CA\t$f_wCpG_CT\t$f_wCpG_CG\t$f_wCpG_GA\t$f_wCpG_GT\t$f_wCpG_GC\t$f_NumDeletion\t$f_NumInsertion\t$f_DeletionBP\t$f_InsertionBP\t$h_wCpG_AT\t$h_wCpG_AG\t$h_wCpG_AC\t$h_wCpG_TA\t$h_wCpG_TC\t$h_wCpG_TG\t$h_wCpG_CA\t$h_wCpG_CT\t$h_wCpG_CG\t$h_wCpG_GA\t$h_wCpG_GT\t$h_wCpG_GC\t$h_NumDeletion\t$h_NumInsertion\t$h_DeletionBP\t$h_InsertionBP\t$f_nCpG_AT\t$f_nCpG_AG\t$f_nCpG_AC\t$f_nCpG_TA\t$f_nCpG_TC\t$f_nCpG_TG\t$f_nCpG_CA\t$f_nCpG_CT\t$f_nCpG_CG\t$f_nCpG_GA\t$f_nCpG_GT\t$f_nCpG_GC\t$h_nCpG_AT\t$h_nCpG_AG\t$h_nCpG_AC\t$h_nCpG_TA\t$h_nCpG_TC\t$h_nCpG_TG\t$h_nCpG_CA\t$h_nCpG_CT\t$h_nCpG_CG\t$h_nCpG_GA\t$h_nCpG_GT\t$h_nCpG_GC\t$wCPG_CombinedAT\t$wCPG_CombinedGC\t$nCpG_CombinedAT\t$nCpG_CombinedGC\t$FALCO_CLASS_MATCH\t$FalcoCLASS_MATCH_NUMBER\t$FalcoCONSENSUS_CLASS_MATCH\t$ConFalcoSex_Chrom_Type\t$ROH" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/revisedSummary/"$s"_ultChromAnnotatedWindowsSummary.txt
}

echo -e "s\tp\tREF_SCAFFOLD\tCOUNT\tA\tC\tG\tT\tCpG\ttGC\ttBASES\tpGC\tuDelSV\tuInsSV\tuInvSV\tuSeg_DelSV\tuSeg_InsSV\tBPuDelSV\tBPuInsSV\tBPuInvSV\tBPuSeg_DelSV\tBPuSeg_InsSV\trDNA\trLTR\trLINE\trSINE\trUnknown\tBPrDNA\tBPrLTR\tBPrLINE\tBPrSINE\tBPrUnknown\tGENES\tSCAFF_SIZE\tSIZE_CAT\tTHRUSH_CHROM_MATCH\tKAKAPO_CHROM_MATCH\tSERIEMA_CHROM_MATCH\tTHRUSH_CLASS_MATCH\tKAKAPO_CLASS_MATCH\tSERIEMA_CLASS_MATCH\tOVERALL_CLASS_MATCH\tCLASS_MATCH_NUMBER\tCONSENSUS_CLASS_MATCH\tf_wCpG_AT\tf_wCpG_AG\tf_wCpG_AC\tf_wCpG_TA\tf_wCpG_TC\tf_wCpG_TG\tf_wCpG_CA\tf_wCpG_CT\tf_wCpG_CG\tf_wCpG_GA\tf_wCpG_GT\tf_wCpG_GC\tf_NumDeletion\tf_NumInsertion\tf_DeletionBP\tf_InsertionBP\th_wCpG_AT\th_wCpG_AG\th_wCpG_AC\th_wCpG_TA\th_wCpG_TC\th_wCpG_TG\th_wCpG_CA\th_wCpG_CT\th_wCpG_CG\th_wCpG_GA\th_wCpG_GT\th_wCpG_GC\th_NumDeletion\th_NumInsertion\th_DeletionBP\th_InsertionBP\tf_nCpG_AT\tf_nCpG_AG\tf_nCpG_AC\tf_nCpG_TA\tf_nCpG_TC\tf_nCpG_TG\tf_nCpG_CA\tf_nCpG_CT\tf_nCpG_CG\tf_nCpG_GA\tf_nCpG_GT\tf_nCpG_GC\th_nCpG_AT\th_nCpG_AG\th_nCpG_AC\th_nCpG_TA\th_nCpG_TC\th_nCpG_TG\th_nCpG_CA\th_nCpG_CT\th_nCpG_CG\th_nCpG_GA\th_nCpG_GT\th_nCpG_GC\twCPG_CombinedAT\twCPG_CombinedGC\tnCpG_CombinedAT\tnCpG_CombinedGC\tFALCO_CLASS_MATCH\tFalcoCLASS_MATCH_NUMBER\tFalcoCONSENSUS_CLASS_MATCH\tConFalcoSex_Chrom_Type\tROH" > /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/revisedSummary/"$s"_ultChromAnnotatedWindowsSummary.txt
COUNT=0
WINDOW_POSITIONS=`cat /scratch/jw5478/FalconProject/Projects/RefGenome8/Synthesis/Windows/Beds/"$s"_16KBph2.haploid.win100KB.bed | perl -pe 's#\b\s+\b#__#g'`
for p in $WINDOW_POSITIONS
do
COUNT=`expr $COUNT + 1`
process_window &
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 28 ]
do
sleep 3s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
wait

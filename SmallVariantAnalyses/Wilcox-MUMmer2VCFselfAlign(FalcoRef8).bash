#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --mem=116GB
#SBATCH -p serial
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=28
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

s=$1

module purge
module load perl
module load gencore/1
module load gencore_variant_detection/1.0

DATE=`date '+ %Y%m%d'`

echo "##fileformat=VCFv4.2" > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".vcf
echo "##fileDate=$DATE" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".vcf
echo "##reference=/scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".vcf

SCAFFOLDS=`grep -P "^[0-9]+" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/MUMmerOutputs/"$s"-v-"$s"100KBhap.snps | awk '{print $11}' | sort -k1,1V | uniq`
cp /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/"$s"_16KBph2.100KB.haploid.fasta
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/"$s"_16KBph2.100KB.haploid.fasta

for n in $SCAFFOLDS
do
LENGTH=`grep -P "^$n" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/"$s"_16KBph2.100KB.haploid.fasta.fai | awk '{print $2}'`
echo "##contig=<ID=$n>" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".vcf
sed -i "s%##contig=<ID=$n>%##contig=<ID=$n,length=$LENGTH>%" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".vcf
done

echo "##INFO=<ID=SNV,Number=0,Type=Flag,Description=\"Indicates that the variant is a SNV.\"" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".vcf
echo "##INFO=<ID=INSERTION,Number=0,Type=Flag,Description=\"Indicates that the variant is an INSERTION releative to refernce genome.\"" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".vcf
echo "##INFO=<ID=DELETION,Number=0,Type=Flag,Description=\"Indicates that the variant is a DELETION relative to the query sequence.\"" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".vcf
echo "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".vcf

#Arrange MUMmer SNVS and Indels output into headerless pseudo-VCF Format
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/TEMP/"$s"/
grep -P "^[0-9]+" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/MUMmerOutputs/"$s"-v-"$s"100KBhap.snps | awk '$7==1 && $8==0{print $0}' | awk 'BEGIN{ OFS = "\t" }{if($2!="." && $3!=".")print $11,$1,".",$2,$3,".",".","SNV",$12,$4,$10;else if($2==".") print $11,$1,".",$2,$3,".",".","INSERTION",$12,$4,$10;else if($3==".") print $11,$1,".",$2,$3,".",".","DELETION",$12,$4,$10}'| awk 'BEGIN { FS=OFS="\t"}$8=="DELETION" && $11==-1{$10=$10+1}$1==p1_1 && $2==p2_1 && $8=="DELETION" && $9==p9_1 && p7_1=="SNV"{next}$1==p1_1 && $2==p2_1 && $8=="DELETION" && $9==p9_1{COUNT_1=COUNT_1+1; POSITION_1=p2_1-COUNT_1; prev_1=$1 OFS POSITION_1 OFS $3 OFS PREV_DEL_1 $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11; PREV_DEL_1=PREV_DEL_1 $4; p2_1=p2_1+1; next}$1==p1_2 && $2==p2_2 && $8=="DELETION" && $9==p9_2 && p7_2=="SNV"{next}$1==p1_2 && $2==p2_2 && $8=="DELETION" && $9==p9_2{COUNT_2=COUNT_2+1; POSITION_2=p2_2-COUNT_2; prev_2=$1 OFS POSITION_2 OFS $3 OFS PREV_DEL_2 $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11; PREV_DEL_2=PREV_DEL_2 $4;p2_2=p2_2+1; next}prev_1!="" && $9 ~ /_1$/{print prev_1}$9 ~ /_1$/{prev_1=$0;p1_1=$1;p2_1=$2+1;p9_1=$9;PREV_DEL_1=$4;p7_1=$7;COUNT_1=0}prev_2!="" && $9 ~ /_2$/{print prev_2}$9 ~ /_2$/{prev_2=$0;p1_2=$1;p2_2=$2+1;p9_2=$9;PREV_DEL_2=$4;p7_2=$7;COUNT_2=0}END{if(prev_1!=""){print prev_1}}END{if(prev_2!=""){print prev_2}}' | awk 'BEGIN { FS=OFS="\t"}$1==p1_1 && $2==p2_1 && $8=="INSERTION" && $9==p9_1 && $10==p10_1 && p7_1=="SNV"{next}$1==p1_1 && $2==p2_1 && $8=="INSERTION" && $9==p9_1 && $10==p10_1 && $11==1{COUNT_1=COUNT_1+1; POSITION_1=p10_1-COUNT_1; prev_1=$1 OFS $2 OFS $3 OFS $4 OFS PREV_INSERT_1 $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS POSITION_1 OFS $11; PREV_INSERT_1=PREV_INSERT_1 $5; p2_1=p2_1; p10_1=p10_1+1; next}$1==p1_1 && $2==p2_1 && $8=="INSERTION" && $9==p9_1 && $10==p10_1 && $11==-1{POSITION_1=p10_1; prev_1=$1 OFS $2 OFS $3 OFS $4 OFS $5 PREV_INSERT_1 OFS $6 OFS $7 OFS $8 OFS $9 OFS POSITION_1 OFS $11;PREV_INSERT_1=$5 PREV_INSERT_1;p2_1=p2_1;p10_1=p10_1+1;next}$1==p1_2 && $2==p2_2 && $8=="INSERTION" && $9==p9_2 && $10==p10_2 && p7_2=="SNV"{next}$1==p1_2 && $2==p2_2 && $8=="INSERTION" && $9==p9_2 && $10==p10_2 && $11==1{COUNT_2=COUNT_2+1; POSITION_2=p10_2-COUNT_2; prev_2=$1 OFS $2 OFS $3 OFS $4 OFS PREV_INSERT_2 $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS POSITION_2 OFS $11; PREV_INSERT_2=PREV_INSERT_2 $5; p2_2=p2_2; p10_2=p10_2+1; next}$1==p1_2 && $2==p2_2 && $8=="INSERTION" && $9==p9_2 && $10==p10_2 && $11==-1{POSITION_2=p10_2; prev_2=$1 OFS $2 OFS $3 OFS $4 OFS $5 PREV_INSERT_2 OFS $6 OFS $7 OFS $8 OFS $9 OFS POSITION_2 OFS $11;PREV_INSERT_2=$5 PREV_INSERT_2;p2_2=p2_2;p10_2=p10_2+1;next}prev_1!="" && $9 ~ /_1$/{print prev_1}$9 ~ /_1$/{prev_1=$0;p1_1=$1;p2_1=$2;p9_1=$9;PREV_INSERT_1=$5;p7_1=$7;p10_1=$10+1;COUNT_1=0}prev_2!="" && $9 ~ /_2$/{print prev_2}$9 ~ /_2$/{prev_2=$0;p1_2=$1;p2_2=$2;p9_2=$9;PREV_INSERT_2=$5;p7_2=$7;p10_2=$10+1;COUNT_2=0}END{if(prev_1!=""){print prev_1}}END{if(prev_2!=""){print prev_2}}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/TEMP/"$s"/"$s"-v-"$s".pseudo_vcf


#Adding in Context Information, i.e. removing .'s
cp /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_16KBph2.diploid.fasta /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/"$s"_16KBph2.diploid.fasta
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/"$s"_16KBph2.diploid.fasta

#Subsample into 100000 fragments
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/TEMP/"$s"/Subsampled/
pVCF_LINES=`wc /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/TEMP/"$s"/"$s"-v-"$s".pseudo_vcf | awk '{print $1}'`
pVCF_SUBSAMPLE=`expr $pVCF_LINES / 100000 + 1`
ITERATIONS=`seq $pVCF_SUBSAMPLE`

for i in $ITERATIONS
do
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 28 ]
do
sleep 3s
JOBS=`jobs | wc | awk '{print $1}'`
done

awk -v ITERATION=$i 'BEGIN{SUBSAMPLE_END=ITERATION*100000; SUBSAMPLE_START=SUBSAMPLE_END-100000; OFS=FS="\t"}NR>SUBSAMPLE_START && NR<=SUBSAMPLE_END{print $0}NR>SUBSAMPLE_END{exit}' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/TEMP/"$s"/"$s"-v-"$s".pseudo_vcf > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/TEMP/"$s"/Subsampled/"$s"-v-"$s".sub"$i".pseudo_vcf &
done

wait

function del_context {
VCF_LINES=`wc $1 | awk '{print $1}'`
DELETION_ID=`grep -v "#" $1 | grep "DELETION" | perl -pe 's#[0-9]+_[1-2]+\s+[0-9]+\s+\.\s+[A-Z]+\s+\.\s+\.\s+\.\s+DELETION\s+([0-9]+_[1-2]+)\s+([0-9]+)\s+(-*[1])#$1_$2_$3#'`
for d in $DELETION_ID;
do
DELETION_SCAFFOLD=`echo $d | perl -pe 's#([0-9]+_[1-2]+)_[0-9]+_-*[1]#$1#'`
DELETION_POS=`echo $d | perl -pe 's#[0-9]+_[1-2]+_([0-9]+)_-*[1]#$1#'`
ORIENTATION=`echo $d | perl -pe 's#[0-9]+_[1-2]+_[0-9]+_(-*[1])#$1#'`
DELETION_AltBASE=`samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/"$s"_16KBph2.diploid.fasta $DELETION_SCAFFOLD:$DELETION_POS-$DELETION_POS | grep -v ">"`
if [[ $ORIENTATION == -1 ]]
then
DELETION_AltBASE=`echo $DELETION_AltBASE | tr ACGTYRWSKMDVHBXNacgtyrwskmdvhbxn TGCARYWSMKHBDVXNtgcarywsmkhbdvxn`
fi
REF_SCAFFOLD=`awk -v DEL_SCAFF=$DELETION_SCAFFOLD -v DEL_POS=$DELETION_POS -v DEL_ALT=$DELETION_AltBASE '$9==DEL_SCAFF && $10==DEL_POS{print $1; exit}' $1`
REF_POS=`awk -v DEL_SCAFF=$DELETION_SCAFFOLD -v DEL_POS=$DELETION_POS -v DEL_ALT=$DELETION_AltBASE '$9==DEL_SCAFF && $10==DEL_POS{print $2; exit}' $1`
REF_prevPOS=`expr $REF_POS - 1`
REF_PREV=`samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/"$s"_16KBph2.100KB.haploid.fasta $REF_SCAFFOLD:$REF_prevPOS-$REF_prevPOS | grep -v ">"`
/home/jw5478/Executables/Programs/gawk-5.1.0/bin/awk -i inplace -v DEL_SCAFF=$DELETION_SCAFFOLD -v DEL_POS=$DELETION_POS -v DEL_ALT=$DELETION_AltBASE -v REF_SCAFF=$REF_SCAFFOLD -v REF_POS=$REF_POS -v REF_PREV=$REF_PREV -v VCF_LINES=VCF_LINES 'BEGIN{OFS="\t"; INDEL_LINES=VCF_LINES}NR>INDEL_LINES{print; next}$1==REF_SCAFF && $2==REF_POS && $9==DEL_SCAFF && $10==DEL_POS{$2=$2-1; $4=REF_PREV$4; $5=DEL_ALT; print; INDEL_LINES=NR; next}{print $0}' $1
done
}

for i in $ITERATIONS
do
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 28 ]
do
sleep 3s
JOBS=`jobs | wc | awk '{print $1}'`
done

del_context /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/TEMP/"$s"/Subsampled/"$s"-v-"$s".sub"$i".pseudo_vcf &
done

wait

function ins_context {
INSERT_ID=`grep -v "#" $1 | grep "INSERTION" | perl -pe 's#[0-9]+_[1-2]+\s+[0-9]+\s+\.\s+\.+\s+[A-Z]+\s+\.\s+\.\s+INSERTION\s+([0-9]+_[1-2]+)\s+([0-9]+)\s+(-*[1])#$1_$2_$3#'`
for d in $INSERT_ID;
do
VCF_LINES=`wc $1 | awk '{print $1}'`
INSERT_SCAFFOLD=`echo $d | perl -pe 's#([0-9]+_[1-2]+)_[0-9]+_-*[1]#$1#'`
INSERT_POS=`echo $d | perl -pe 's#[0-9]+_[1-2]+_([0-9]+)_-*[1]#$1#'`
ORIENTATION=`echo $d | perl -pe 's#[0-9]+_[1-2]+_[0-9]+_(-*[1])#$1#'`
if [[ $ORIENTATION == 1 ]]
then
INSERT_prevPOS=`expr $INSERT_POS - 1`
INSERT_PREV=`samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/"$s"_16KBph2.diploid.fasta $INSERT_SCAFFOLD:$INSERT_prevPOS-$INSERT_prevPOS | grep -v ">"`
fi
if [[ $ORIENTATION == -1 ]]
then
INSERT_prevPOS=`expr $INSERT_POS + 1`
INSERT_PREV=`samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/"$s"_16KBph2.diploid.fasta $INSERT_SCAFFOLD:$INSERT_prevPOS-$INSERT_prevPOS | grep -v ">" | tr ACGTYRWSKMDVHBXNacgtyrwskmdvhbxn TGCARYWSMKHBDVXNtgcarywsmkhbdvxn`
fi
REF_SCAFFOLD=`awk -v INS_SCAFF=$INSERT_SCAFFOLD -v INS_POS=$INSERT_POS -v INS_PREV=$INSERT_PREV '$9==INS_SCAFF && $10==INS_POS{print $1; exit}' $1`
REF_POS=`awk -v INS_SCAFF=$INSERT_SCAFFOLD -v INS_POS=$INSERT_POS -v INS_PREV=$INSERT_PREV '$9==INS_SCAFF && $10==INS_POS{print $2; exit}' $1`
REF_BASE=`samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/"$s"_16KBph2.100KB.haploid.fasta $REF_SCAFFOLD:$REF_POS-$REF_POS | grep -v ">"`
/home/jw5478/Executables/Programs/gawk-5.1.0/bin/awk -i inplace -v INS_SCAFF=$INSERT_SCAFFOLD -v INS_POS=$INSERT_POS -v INS_PREV=$INSERT_PREV -v REF_SCAFF=$REF_SCAFFOLD -v REF_POS=$REF_POS -v REF_BASE=$REF_BASE -v VCF_LINES=VCF_LINES 'BEGIN{OFS="\t"; INDEL_LINES=VCF_LINES}NR>INDEL_LINES{print; next}$9==INS_SCAFF && $10==INS_POS && $1==REF_SCAFF && $2==REF_POS && $11==1{$4=REF_BASE; $5=INS_PREV$5; $10=$10-1; print; INDEL_LINES=NR; next}$9==INS_SCAFF && $10==INS_POS && $1==REF_SCAFF && $2==REF_POS && $11==-1{$4=REF_BASE; $5=INS_PREV$5; $10=$10+1; print; INDEL_LINES=NR; next}{print $0}' $1
done
}

for i in $ITERATIONS
do
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 28 ]
do
sleep 3s
JOBS=`jobs | wc | awk '{print $1}'`
done

ins_context /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/TEMP/"$s"/Subsampled/"$s"-v-"$s".sub"$i".pseudo_vcf &
done

wait

rm /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".context.pseudo_vcf
rm /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".context.sort.pseudo_vcf
for i in $ITERATIONS
do
cat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/TEMP/"$s"/Subsampled/"$s"-v-"$s".sub"$i".pseudo_vcf >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".context.pseudo_vcf
done
cat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".context.pseudo_vcf | sort -k1,1V -k2,2n > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".context.sort.pseudo_vcf
awk 'BEGIN{ FS=OFS="\t"}{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8}' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".context.sort.pseudo_vcf >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".vcf
bcftools norm -f /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/"$s"_16KBph2.100KB.haploid.fasta -c ws /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".vcf > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".norm.vcf

bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".vcf &
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".norm.vcf &
wait

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

DATE=`date '+ %Y%m%d'`

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques
bcftools isec --nfiles=1 -c none -p /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/*-v-Kestrel100KBhap.dip.norm.vcf.gz

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts

echo  "$s\tINSERT\tDEL\tAT\tAC\tAG\tTA\tTC\tTG\tCA\tCT\tCG\tGA\tGT\tGC" | perl -pe "s#\\\t#\t#g" > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/UniqueKestrelAlignVariantCount.txt
SEQs=`seq 0 7`
for s in $SEQs
do
#'Insertion'
INSERT=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk '$8=="INSERTION"{print $0}' | wc | awk '{print $1}'`
#'Deletion'
DEL=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk '$8=="DELETION"{print $0}' | wc | awk '{print $1}'`

#'A>T'
AT=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk '$4=="A" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'A>C'
AC=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk '$4=="A" && $5=="C"{print $0}' | wc | awk '{print $1}'`
#'A>G'
AG=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk '$4=="A" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'T>A'
TA=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk '$4=="T" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'T>C'
TC=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk '$4=="T" && $5=="C"{print $0}' | wc | awk '{print $1}'`
#'T>G'
TG=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk '$4=="T" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'C>A'
CA=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk '$4=="C" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'C>T'
CT=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk '$4=="C" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'C>G'
CG=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk '$4=="C" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'G>A'
GA=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk '$4=="G" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'G>T'
GT=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk '$4=="G" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'G>C'
GC=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk '$4=="G" && $5=="C"{print $0}' | wc | awk '{print $1}'`

s=`grep 000"$s".vcf /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/README.txt | perl -pe 's#.*/VCFs/Diploid/##' | perl -pe 's#-v-Kestrel100KBhap.dip.norm.vcf.gz##'`

echo "$s\t$INSERT\t$DEL\t$AT\t$AC\t$AG\t$TA\t$TC\t$TG\t$CA\t$CT\t$CG\t$GA\t$GT\t$GC" | perl -pe "s#\\\t#\t#g" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/UniqueKestrelAlignVariantCount.txt
done

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs


echo "##fileformat=VCFv4.2" > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs/Kestrel_16KBph2.100KB.haploid.gvcf
echo "##fileDate=$DATE" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs/Kestrel_16KBph2.100KB.haploid.gvcf
echo "##reference=/scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/Kestrel_16KBph2.100KB.haploid.fasta" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs/Kestrel_16KBph2.100KB.haploid.gvcf

SCAFFOLDS=`awk '{print $1}' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Kestrel_16KBph2.100KB.haploid.fasta.fai | sort -k1,1V | uniq`

for n in $SCAFFOLDS
do
LENGTH=`grep -P "^$n" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Kestrel_16KBph2.100KB.haploid.fasta.fai | awk '{print $2}'`
echo "##contig=<ID=$n>" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs/Kestrel_16KBph2.100KB.haploid.gvcf
sed -i "s%##contig=<ID=$n>%##contig=<ID=$n,length=$LENGTH>%" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs/Kestrel_16KBph2.100KB.haploid.gvcf
done

echo "##INFO=<ID=SNV,Number=0,Type=Flag,Description=\"Indicates that the variant is a SNV.\"" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs/Kestrel_16KBph2.100KB.haploid.gvcf
echo "##INFO=<ID=INSERTION,Number=0,Type=Flag,Description=\"Indicates that the variant is an INSERTION releative to refernce genome.\"" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs/Kestrel_16KBph2.100KB.haploid.gvcf
echo "##INFO=<ID=DELETION,Number=0,Type=Flag,Description=\"Indicates that the variant is a DELETION relative to the query sequence.\"" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs/Kestrel_16KBph2.100KB.haploid.gvcf
echo "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs/Kestrel_16KBph2.100KB.haploid.gvcf

perl -pe 's#(>[0-9]+_[1-2]).*#$1#' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Kestrel_16KBph2.100KB.haploid.fasta | perl -pe 's#([a-zA-Z]+)\n#$1#' | perl -pe 's#([a-zA-Z])#$1\n#g' | awk 'BEGIN{OFS=FS="\t"}/>/{Length=length($0); Scaff=substr($0,2,Length); Count=1; next}{print Scaff OFS  Count OFS "." OFS $1 OFS $ 1 OFS "." OFS "." OFS "REF"; Count=Count+1}' >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs/Kestrel_16KBph2.100KB.haploid.gvcf

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/vReference

bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs/Kestrel_16KBph2.100KB.haploid.gvcf
bcftools index /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs/Kestrel_16KBph2.100KB.haploid.gvcf.gz
bcftools isec --nfiles=1 -c all -w 1 -p /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/vReference /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/gVCFs/Kestrel_16KBph2.100KB.haploid.gvcf.gz /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/*-v-Kestrel100KBhap.dip.norm.vcf.gz

rm /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/ConservedSitesBaseComp.txt
echo -e Conserved_Bases'\t'A'\t'T'\t'G'\t'C'\t'CombinedAT'\t'CombinedGC > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/HaploidConservedSitesBaseComp.txt
echo -e Haploid_Count'\t'$A'\t'$T'\t'$G'\t'$C'\t'$CombinedAT'\t'$CombinedGC >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/HaploidConservedSitesBaseComp.txt

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC
rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC

#With no GC
##Extract any sites that are CpG in any of the genomes into a single vcf
for s in A D B FpBSx77 C E F Kest2
do
zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf.gz | grep "#" > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/"$s"_16KBph2.Kest100KBhap.dip_con.GC.vcf
zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf.gz | grep -v "#" | awk 'BEGIN{FS=OFS="\t"}$1==scaff && $2==pos + 1 && ($4==check_base || $5==check_base){print prev "\n" $0}$4=="C" || $5=="C"{scaff=$1; pos=$2; prev=$0; check_base="G"}$4=="G" || $5=="G"{scaff=$1; pos=$2; prev=$0; check_base="C"}' >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/"$s"_16KBph2.Kest100KBhap.dip_con.GC.vcf
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/"$s"_16KBph2.Kest100KBhap.dip_con.GC.vcf
done

zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/Kest2_16KBph2.Kest100KBhap.dip_con.gvcf.gz | grep "#" > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/AllRef-GC_16KBph2.Kest100KBhap.dip_con.GC.vcf
zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/*_16KBph2.Kest100KBhap.dip_con.GC.vcf.gz | grep -v "#" | sort -k1,1V -k2,2n | uniq >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/AllRef-GC_16KBph2.Kest100KBhap.dip_con.GC.vcf


rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/vReference
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/AllRef-GC_16KBph2.Kest100KBhap.dip_con.GC.vcf
bcftools index /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/AllRef-GC_16KBph2.Kest100KBhap.dip_con.GC.vcf.gz
cp /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/vReference/0000.vcf /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/ConservedSites.vcf
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/ConservedSites.vcf
bcftools index /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/ConservedSites.vcf.gz
bcftools isec --nfiles=1 -c all -w 1 -p /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/vReference /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/ConservedSites.vcf.gz /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/AllRef-GC_16KBph2.Kest100KBhap.dip_con.GC.vcf.gz


rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs
rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/Rename
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/Rename
cp /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/*.vcf /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/Rename/

for n in $(ls /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/Rename)
do
NAME=`grep $n /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/README.txt | perl -pe 's#.+/(.+)-v-Kestrel100KBhap.dip.norm.vcf.gz#$1#g'`
mv /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/Rename/$n  /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/Rename/"$NAME".unique.vcf
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/Rename/"$NAME".unique.vcf
bcftools index /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/Rename/"$NAME".unique.vcf.gz
done

for s in A B C D E FpBSx77 F Kestrel
do
bcftools isec --nfiles=1 -c all -w 1 -p /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/Rename/"$s".unique.vcf.gz /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/AllRef-GC_16KBph2.Kest100KBhap.dip_con.GC.vcf.gz
done

#Summary
##Conserved Sites Summary
A=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/vReference/0000.vcf | awk 'BEGIN{Count=0}$4=="A"{Count=Count+1}END{print Count}'`
T=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/vReference/0000.vcf | awk 'BEGIN{Count=0}$4=="T"{Count=Count+1}END{print Count}'`
G=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/vReference/0000.vcf | awk 'BEGIN{Count=0}$4=="G"{Count=Count+1}END{print Count}'`
C=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/vReference/0000.vcf | awk 'BEGIN{Count=0}$4=="C"{Count=Count+1}END{print Count}'`
CombinedAT=`expr $A + $T`
CombinedGC=`expr $G + $C`

rm /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/NoGC_ConservedSitesBaseComp.txt
echo -e Conserved_Bases'\t'A'\t'T'\t'G'\t'C'\t'CombinedAT'\t'CombinedGC > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/NoGC_ConservedSitesBaseComp.txt
echo -e Haploid_Count'\t'$A'\t'$T'\t'$G'\t'$C'\t'$CombinedAT'\t'$CombinedGC >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/NoGC_ConservedSitesBaseComp.txt

#Summary of Unique Variants
rm /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/NoGC_UniqueKestrelAlignVariantCount.txt
echo  "$s\tINSERT\tDEL\tAT\tAC\tAG\tTA\tTC\tTG\tCA\tCT\tCG\tGA\tGT\tGC" | perl -pe "s#\\\t#\t#g" > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/NoGC_UniqueKestrelAlignVariantCount.txt
for s in A B C D E FpBSx77 F Kestrel
do
#'Insertion'
INSERT=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s"/0000.vcf | awk '$8=="INSERTION"{print $0}' | wc | awk '{print $1}'`
#'Deletion'
DEL=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s"/0000.vcf | awk '$8=="DELETION"{print $0}' | wc | awk '{print $1}'`

#'A>T'
AT=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s"/0000.vcf | awk '$4=="A" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'A>C'
AC=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s"/0000.vcf | awk '$4=="A" && $5=="C"{print $0}' | wc | awk '{print $1}'`
#'A>G'
AG=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s"/0000.vcf | awk '$4=="A" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'T>A'
TA=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s"/0000.vcf | awk '$4=="T" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'T>C'
TC=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s"/0000.vcf | awk '$4=="T" && $5=="C"{print $0}' | wc | awk '{print $1}'`
#'T>G'
TG=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s"/0000.vcf | awk '$4=="T" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'C>A'
CA=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s"/0000.vcf | awk '$4=="C" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'C>T'
CT=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s"/0000.vcf | awk '$4=="C" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'C>G'
CG=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s"/0000.vcf | awk '$4=="C" && $5=="G"{print $0}' | wc | awk '{print $1}'`

#'G>A'
GA=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s"/0000.vcf | awk '$4=="G" && $5=="A"{print $0}' | wc | awk '{print $1}'`
#'G>T'
GT=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s"/0000.vcf | awk '$4=="G" && $5=="T"{print $0}' | wc | awk '{print $1}'`
#'G>C'
GC=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/noGC/SNVs/"$s"/0000.vcf | awk '$4=="G" && $5=="C"{print $0}' | wc | awk '{print $1}'`

echo "$s\t$INSERT\t$DEL\t$AT\t$AC\t$AG\t$TA\t$TC\t$TG\t$CA\t$CT\t$CG\t$GA\t$GT\t$GC" | perl -pe "s#\\\t#\t#g" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/NoGC_UniqueKestrelAlignVariantCount.txt
done

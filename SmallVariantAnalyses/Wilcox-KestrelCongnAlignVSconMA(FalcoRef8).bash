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

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/Phased
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Alignments

cp /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/vReference/0000.vcf /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/ConservedSites.vcf
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/ConservedSites.vcf
bcftools index /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/ConservedSites.vcf.gz

function consensus_vcf {
zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/"$s"-v-Kestrel100KBhap.dip.norm.vcf.gz | grep "#" > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/"$s"-v-Kestrel100KBhap.dip.norm.snv.vcf
zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/"$s"-v-Kestrel100KBhap.dip.norm.vcf.gz | grep -v "#" | grep "SNV" | awk 'length($4)==1 && length($5)==1{print $0}' >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/"$s"-v-Kestrel100KBhap.dip.norm.snv.vcf
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/"$s"-v-Kestrel100KBhap.dip.norm.snv.vcf
bcftools index /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/"$s"-v-Kestrel100KBhap.dip.norm.snv.vcf.gz

zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/"$s"-v-Kestrel100KBhap.dip.norm.snv.vcf.gz | grep "#" > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/Phased/"$s"-v-Kestrel100KBhap.dip.norm.snv.phased.vcf
line_count=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/"$s"-v-Kestrel100KBhap.dip.norm.snv.vcf.gz | grep -v -c "#"`
zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/"$s"-v-Kestrel100KBhap.dip.norm.snv.vcf.gz | grep -v "#" | awk -v final=$line_count 'BEGIN{FS=OFS="\t"; new="yes"}new=="yes" && NR==final{print $1 OFS $2 OFS $3 OFS $4 OFS $4 "|" $5 OFS $6 OFS $7 OFS $8; exit}($1!=chrom || $2!=pos) && new=="yes"{chrom=$1; pos=$2; alt=$5; snv=chrom OFS pos OFS "." OFS $4 OFS $4 "|" $5 OFS $6 OFS $7 OFS $8; new="no"; next}$1==chrom && $2==pos && new=="no"{print chrom OFS pos OFS "." OFS $4 OFS alt "|" $5 OFS $6 OFS $7 OFS $8; new="yes"; next}($1!=chrom || $2!=pos) && new=="no"{print snv; new="yes"}' | perl -pe 's#A\|C#M#' | perl -pe 's#C\|A#M#' | perl -pe 's#A\|G#R#' | perl -pe 's#G\|A#R#' | perl -pe 's#A\|T#W#' | perl -pe 's#T\|A#W#' | perl -pe 's#C\|G#S#' | perl -pe 's#G\|C#S#' | perl -pe 's#C\|T#Y#' | perl -pe 's#T\|C#Y#' | perl -pe 's#G\|T#K#' | perl -pe 's#T\|G#K#' | perl -pe 's#\|[a-zA-Z]##' >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/Phased/"$s"-v-Kestrel100KBhap.dip.norm.snv.phased.vcf
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/Phased/"$s"-v-Kestrel100KBhap.dip.norm.snv.phased.vcf
bcftools index /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/Phased/"$s"-v-Kestrel100KBhap.dip.norm.snv.phased.vcf.gz
cat /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/Kestrel_16KBph2.100KB.haploid.fasta | bcftools consensus /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/Phased/"$s"-v-Kestrel100KBhap.dip.norm.snv.phased.vcf.gz > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/"$s"_Kest100KBhap.dip_con.fasta
}

for s in A D B FpBSx77 C E F
do
consensus_vcf &
done

zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/Kestrel-v-Kestrel100KBhap.dip.norm.vcf.gz  | grep "#" > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/Kestrel-v-Kestrel100KBhap.dip.norm.snv.vcf
zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/Kestrel-v-Kestrel100KBhap.dip.norm.vcf.gz  | grep -v "#" | grep "SNV" | awk 'length($4)==1 && length($5)==1{print $0}' | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $3 OFS $4 OFS $4"|"$5 OFS $6 OFS $7 OFS $8}' | perl -pe 's#A\|C#M#' | perl -pe 's#C\|A#M#' | perl -pe 's#A\|G#R#' | perl -pe 's#G\|A#R#' | perl -pe 's#A\|T#W#' | perl -pe 's#T\|A#W#' | perl -pe 's#C\|G#S#' | perl -pe 's#G\|C#S#' | perl -pe 's#C\|T#Y#' | perl -pe 's#T\|C#Y#' | perl -pe 's#G\|T#K#' | perl -pe 's#T\|G#K#' >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/Kestrel-v-Kestrel100KBhap.dip.norm.snv.vcf
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/Kestrel-v-Kestrel100KBhap.dip.norm.snv.vcf
bcftools index /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/Kestrel-v-Kestrel100KBhap.dip.norm.snv.vcf.gz
cat /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/Kestrel_16KBph2.100KB.haploid.fasta | bcftools consensus /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/VCFs/Kestrel-v-Kestrel100KBhap.dip.norm.snv.vcf.gz > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/Kest2_Kest100KBhap.dip_con.fasta
wait

#Make gVCFs of each consensus
DATE=`date '+ %Y%m%d'`

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs

function make_gvcf {
echo "##fileformat=VCFv4.2" > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf
echo "##fileDate=$DATE" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf
echo "##reference=/scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf

SCAFFOLDS=`awk '{print $1}' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Kestrel_16KBph2.100KB.haploid.fasta.fai | sort -k1,1V | uniq`

for n in $SCAFFOLDS
do
LENGTH=`grep -P "^$n" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Kestrel_16KBph2.100KB.haploid.fasta.fai | awk '{print $2}'`
echo "##contig=<ID=$n>" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf
sed -i "s%##contig=<ID=$n>%##contig=<ID=$n,length=$LENGTH>%" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf
done

echo "##INFO=<ID=SNV,Number=0,Type=Flag,Description=\"Indicates that the variant is a SNV.\"" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf
echo "##INFO=<ID=INSERTION,Number=0,Type=Flag,Description=\"Indicates that the variant is an INSERTION releative to refernce genome.\"" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf
echo "##INFO=<ID=DELETION,Number=0,Type=Flag,Description=\"Indicates that the variant is a DELETION relative to the query sequence.\"" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf
echo "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf

perl -pe 's#(>[0-9]+_[1-2]).*#$1#' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/"$s"_Kest100KBhap.dip_con.fasta | perl -pe 's#([a-zA-Z])#$1\n#g' | awk 'BEGIN{OFS=FS="\t"}/>/{Length=length($0); Scaff=substr($0,2,Length); Count=1; next}$0!=""{print Scaff OFS  Count OFS "." OFS $1 OFS $1 OFS "." OFS "." OFS "REF"; Count=Count+1}' >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf
bcftools index /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf.gz
}
for s in A D B FpBSx77 C E F Kest2
do
make_gvcf &
done
wait

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/VarialbeSites
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/VarialbeSites

for s in A D B FpBSx77 C E F Kest2
do
bcftools isec --nfiles=1 -c none -w 1 -p /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/VarialbeSites/"$s" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/"$s"_16KBph2.Kest100KBhap.dip_con.gvcf.gz  /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/GC/ConservedSites.vcf.gz &
done
wait

function bgzip_it {
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/VarialbeSites/"$s"/0000.vcf 
bcftools index /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/VarialbeSites/"$s"/0000.vcf.gz
}

for s in A D B FpBSx77 C E F Kest2
do
bgzip_it &
done
wait

rm /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Alignments/MA_Kest100KBhap.VarSite_dip_con.vs.fasta
for s in A D B FpBSx77 C E F Kest2
do
echo ">$s.genome" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Alignments/MA_Kest100KBhap.VarSite_dip_con.vs.fasta
zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Genomes/gVCFs/VarialbeSites/"$s"/0000.vcf.gz | grep -v "^#" | awk '{print $5}' | perl -pe -chomp >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Alignments/MA_Kest100KBhap.VarSite_dip_con.vs.fasta
echo -e "\n" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Alignments/MA_Kest100KBhap.VarSite_dip_con.vs.fasta
done

java -cp /home/jw5478/Executables/Programs/readseq.jar run -inform=8 -f=17 -p /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Alignments/MA_Kest100KBhap.VarSite_dip_con.vs.fasta > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Alignments/MA_Kest100KBhap.VarSite_dip_con.vs.nexus
java -cp /home/jw5478/Executables/Programs/readseq.jar run -inform=8 -f=12 -p /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Alignments/MA_Kest100KBhap.VarSite_dip_con.vs.fasta > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/Consensus/Alignments/MA_Kest100KBhap.VarSite_dip_con.vs.phylip

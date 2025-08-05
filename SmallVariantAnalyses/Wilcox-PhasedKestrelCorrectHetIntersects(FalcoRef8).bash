module purge
module load perl
module load gencore/1
module load gencore_variant_detection/1.0

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/SharedSites
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/SharedSNV-INDELs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/SharedSNV-INDELs/SharedSNVs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/SharedSNV-INDELs/SharedIndels

mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/

for s in A D B FpBSx77 C E F
do
zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/"$s"-v-Kestrel100KBhap.dip.vcf.gz | grep "#" > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap1.vcf
awk 'BEGIN{ FS=OFS="\t"}$9~/_1/{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8}' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/TEMP/"$s"/"$s"-v-Kestrel100KBhap.dip.context.sort.pseudo_vcf >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap1.vcf
bcftools norm -f /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/Kestrel_16KBph2.100KB.haploid.fasta -c ws /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap1.vcf > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap1.norm.vcf &

zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/"$s"-v-Kestrel100KBhap.dip.vcf.gz | grep "#" > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap2.vcf
awk 'BEGIN{ FS=OFS="\t"}$9~/_2/{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8}' /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/TEMP/"$s"/"$s"-v-Kestrel100KBhap.dip.context.sort.pseudo_vcf >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap2.vcf
bcftools norm -f /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/Kestrel_16KBph2.100KB.haploid.fasta -c ws /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap2.vcf > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap2.norm.vcf &

wait

bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap1.vcf &
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap1.norm.vcf &

bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap2.vcf &
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap2.norm.vcf &

wait

bcftools index -f /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap1.norm.vcf.gz &
bcftools index -f /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap2.norm.vcf.gz &

wait

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$s"
rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/Kest_Hap"$s"-v-"$s"/

bcftools isec -c none -p /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/Kest_Hap"$s"-v-"$s"/ /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap1.norm.vcf.gz /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Phased/"$s"-v-Kestrel100KBhap.hap2.norm.vcf.gz

bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/Kest_Hap"$s"-v-"$s"/0000.vcf 
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/Kest_Hap"$s"-v-"$s"/0001.vcf

bcftools index /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/Kest_Hap"$s"-v-"$s"/0000.vcf.gz
bcftools index /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/Kest_Hap"$s"-v-"$s"/0001.vcf.gz

bcftools isec --nfiles=1 -c none -p /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$s"/ /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/Kest_Hap"$s"-v-"$s"/0000.vcf.gz /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/Kest_Hap"$s"-v-"$s"/0001.vcf.gz /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/Kestrel-v-Kestrel100KBhap.dip.norm.vcf.gz


bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$s"/0000.vcf &
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$s"/0001.vcf &

wait

bcftools index -f /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$s"/0000.vcf.gz &
bcftools index -f /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$s"/0001.vcf.gz &

wait

bcftools isec -c none -p /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/KestrelDipStrip/"$s" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/"$s"-v-Kestrel100KBhap.dip.norm.vcf.gz /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/Kestrel-v-Kestrel100KBhap.dip.norm.vcf.gz
bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/KestrelDipStrip/"$s"/0000.vcf
bcftools index /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/KestrelDipStrip/"$s"/0000.vcf.gz
done

COMPARATIVE=`echo D B FpBSx77 C E F`
for s in A D B FpBSx77 C E F

do
for c in $COMPARATIVE
do
rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/SharedSNV-INDELs/"$s"-v-"$c"/
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/SharedSNV-INDELs/"$s"-v-"$c"/
rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/SharedSNV-INDELs/"$s"-v-"$c"/KestDipStrip
bcftools isec -c none -p /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/SharedSNV-INDELs/"$s"-v-"$c"/KestDipStrip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/KestrelDipStrip/"$s"/0000.vcf.gz /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/KestrelDipStrip/"$c"/0000.vcf.gz

bgzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/SharedSNV-INDELs/"$s"-v-"$c"/KestDipStrip/0000.vcf
bcftools index /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/SharedSNV-INDELs/"$s"-v-"$c"/KestDipStrip/0000.vcf.gz

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$c"
bcftools isec --nfiles=2 -c none -p /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$c" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/SharedSNV-INDELs/"$s"-v-"$c"/KestDipStrip/0000.vcf.gz /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$s"/0000.vcf.gz /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$s"/0001.vcf.gz 

bcftools isec --nfiles=2 -c none -p /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$c"/"$c"_Intersect /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/SharedSNV-INDELs/"$s"-v-"$c"/KestDipStrip/0000.vcf.gz /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$c"-v-"$c"/0000.vcf.gz /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$c"-v-"$c"/0001.vcf.gz

done
COMPARATIVE=`echo $COMPARATIVE | perl -pe "s#[a-zA-Z0-9]+\s*##"`
done

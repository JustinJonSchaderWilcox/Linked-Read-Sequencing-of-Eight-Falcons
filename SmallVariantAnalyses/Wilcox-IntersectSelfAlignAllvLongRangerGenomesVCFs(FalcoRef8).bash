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

module purge
module load perl
module load gencore/1
module load gencore_variant_detection/1.0

DATE=`date '+ %Y%m%d'`

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs
rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm
rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/Kestrel-SelfVCFs/
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/Kestrel-SelfVCFs/
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp

for s in A B C D E F FpBSx77 Kestrel
do
#Normalize Longranger outputs for comparison
zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".norm.vcf.gz | grep "#" | grep -v "##contig=" | grep -v "#CHROM" | bgzip > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.vcf.gz
zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/LongrangerOut/SelfAlign/16KBph2Haploid/"$s"_16KBph2-HaploidLR/outs/phased_variants.vcf.gz | grep "##contig=<ID=" | bgzip >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.vcf.gz
zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".norm.vcf.gz | grep "#CHROM" | bgzip >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.vcf.gz
zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/LongrangerOut/SelfAlign/16KBph2Haploid/"$s"_16KBph2-HaploidLR/outs/phased_variants.vcf.gz | awk 'BEGIN{ OFS = "\t" }{if($7=="PASS")print $0}' | perl -pe 's#(\s)AC=.*$#$1\.#' | bgzip >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.vcf.gz
bcftools norm -f /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta -c ws /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.vcf.gz | bgzip > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.norm.vcf.gz

#Index VCFs
bcftools index -f /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".norm.vcf.gz
bcftools index -f /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.norm.vcf.gz

mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/"$s"

#Shared SNV-INDELs within each Genome
##Longranger versus MUMmer self align
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/"$s"/LR-v-MUMmerSA
bcftools isec -c none -p /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/"$s"/LR-v-MUMmerSA /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".norm.vcf.gz /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.norm.vcf.gz
done

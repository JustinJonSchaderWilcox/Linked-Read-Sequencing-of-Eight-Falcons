rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/Summary
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/Summary

echo -e Genome'\t'mumHET_SNVs'\t'mumHET_INDELs'\t'lrHET_SNVs'\t'lrHET_INDELs'\t'UmumHET_SNVs'\t'UmumHET_INDELs'\t'UlrHET_SNVs'\t'UlrHET_INDELs'\t'KEST_SNVsHet'\t'KEST_INDELs > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/Summary/MethodsCompVariantSummary.txt

for s in A D B FpBSx77 C E F
do
mumHET_SNVs=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".norm.vcf.gz | grep -c "SNV"`
mumHET_DELs=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".norm.vcf.gz | grep -c "DELETION"`
mumHET_INSERTIONs=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".norm.vcf.gz | grep -c "INSERTION"`
mumHET_INDELs=`expr $mumHET_DELs + $mumHET_INSERTIONs`

lrHET_SNVs=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.norm.vcf.gz | grep -v "#" | awk '{lRef=length($4); lAlt=length($5)}lRef==1 && lAlt==1{print "SNV"}' | grep -c "SNV"`
lrHET_INDELs=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.norm.vcf.gz | grep -v "#" | awk '{lRef=length($4); lAlt=length($5)}lRef!=1 || lAlt!=1{print "InDel"}' | grep -c "InDel"`

UmumHET_SNVs=`grep -c "SNV" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/$s/LR-v-MUMmerSA/0000.vcf`
UmumHET_DELs=`grep -c "DELETION" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/$s/LR-v-MUMmerSA/0000.vcf`
UmumHET_INSERTIONs=`grep -c "INSERTION" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/$s/LR-v-MUMmerSA/0000.vcf`
UmumHET_INDELs=`expr $UmumHET_DELs + $UmumHET_INSERTIONs`

UlrHET_SNVs=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/$s/LR-v-MUMmerSA/0001.vcf | awk '{lRef=length($4); lAlt=length($5)}lRef==1 && lAlt==1{print "SNV"}' | grep -c "SNV"`
UlrHET_INDELs=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/$s/LR-v-MUMmerSA/0001.vcf | awk '{lRef=length($4); lAlt=length($5)}lRef!=1 || lAlt!=1{print "InDel"}' | grep -c "InDel"`

DIV_SNVsHap1=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$s"/0000.vcf.gz | grep -v "#" | grep -c "SNV"`
DIV_SNVsHap2=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$s"/0001.vcf.gz | grep -v "#" | grep -c "SNV"`
KEST_SNVsHet=`expr $DIV_SNVsHap1 + $DIV_SNVsHap2`

DIV_DELsHap1=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$s"/0000.vcf.gz | grep -v "#" | grep -c "DELETION"`
DIV_DELsHap2=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$s"/0001.vcf.gz | grep -v "#" | grep -c "DELETION"`
DIV_DELsHet=`expr $DIV_DELsHap1 + $DIV_DELsHap2`

DIV_INSERTIONsHap1=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$s"/0000.vcf.gz | grep -v "#" | grep -c "INSERTION"`
DIV_INSERTIONsHap2=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Heterozygotes/SNV-INDELs/"$s"-v-"$s"/0001.vcf.gz | grep -v "#" | grep -c "INSERTION"`
DIV_INSERTIONsHet=`expr $DIV_INSERTIONsHap1 + $DIV_INSERTIONsHap2`

KEST_INDELs=`expr $DIV_DELsHet + $DIV_INSERTIONsHet`

echo -e $s'\t'$mumHET_SNVs'\t'$mumHET_INDELs'\t'$lrHET_SNVs'\t'$lrHET_INDELs'\t'$UmumHET_SNVs'\t'$UmumHET_INDELs'\t'$UlrHET_SNVs'\t'$UlrHET_INDELs'\t'$KEST_SNVsHet'\t'$KEST_INDELs >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/Summary/MethodsCompVariantSummary.txt
done

s=Kestrel
mumHET_SNVs=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".norm.vcf.gz | grep -c "SNV"`
mumHET_DELs=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".norm.vcf.gz | grep -c "DELETION"`
mumHET_INSERTIONs=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/VCFs/"$s"-v-"$s".norm.vcf.gz | grep -c "INSERTION"`
mumHET_INDELs=`expr $mumHET_DELs + $mumHET_INSERTIONs`

lrHET_SNVs=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.norm.vcf.gz | grep -v "#" | awk '{lRef=length($4); lAlt=length($5)}lRef==1 && lAlt==1{print "SNV"}' | grep -c "SNV"`
lrHET_INDELs=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/LongrangerNorm/"$s"_LR_variants.norm.vcf.gz | grep -v "#" | awk '{lRef=length($4); lAlt=length($5)}lRef!=1 || lAlt!=1{print "InDel"}' | grep -c "InDel"`

UmumHET_SNVs=`grep -c "SNV" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/$s/LR-v-MUMmerSA/0000.vcf`
UmumHET_DELs=`grep -c "DELETION" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/$s/LR-v-MUMmerSA/0000.vcf`
UmumHET_INSERTIONs=`grep -c "INSERTION" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/$s/LR-v-MUMmerSA/0000.vcf`
UmumHET_INDELs=`expr $UmumHET_DELs + $UmumHET_INSERTIONs`

UlrHET_SNVs=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/$s/LR-v-MUMmerSA/0001.vcf | awk '{lRef=length($4); lAlt=length($5)}lRef==1 && lAlt==1{print "SNV"}' | grep -c "SNV"`
UlrHET_INDELs=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/$s/LR-v-MUMmerSA/0001.vcf | awk '{lRef=length($4); lAlt=length($5)}lRef!=1 || lAlt!=1{print "InDel"}' | grep -c "InDel"`

KEST_SNVsHet=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/Kestrel-v-Kestrel100KBhap.dip.norm.vcf.gz | grep -v "#" | grep -c "SNV"`

DIV_DELsHet=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/Kestrel-v-Kestrel100KBhap.dip.norm.vcf.gz | grep -v "#" | grep -c "DELETION"`
DIV_INSERTIONsHet=`zcat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/VCFs/Diploid/Kestrel-v-Kestrel100KBhap.dip.norm.vcf.gz | grep -v "#" | grep -c "INSERTION"`

KEST_INDELs=`expr $DIV_DELsHet + $DIV_INSERTIONsHet`

echo -e $s'\t'$mumHET_SNVs'\t'$mumHET_INDELs'\t'$lrHET_SNVs'\t'$lrHET_INDELs'\t'$UmumHET_SNVs'\t'$UmumHET_INDELs'\t'$UlrHET_SNVs'\t'$UlrHET_INDELs'\t'$KEST_SNVsHet'\t'$KEST_INDELs >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/SelfAlign/SNV-INDELs/MethodsComp/Summary/MethodsCompVariantSummary.txt

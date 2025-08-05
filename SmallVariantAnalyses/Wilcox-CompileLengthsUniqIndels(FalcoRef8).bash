SEQs=`seq 0 7`
echo  "Genome\tType\tLength" | perl -pe "s#\\\t#\t#g" > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/UniqueKestrelAlignIndelSize.txt
for s in $SEQs
do
g=`grep 000"$s".vcf /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/README.txt | perl -pe 's#.*/VCFs/Diploid/##' | perl -pe 's#-v-Kestrel100KBhap.dip.norm.vcf.gz##'`
grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/000"$s".vcf | awk -v genome="$g" 'BEGIN{OFS=FS="\t"}$8=="INSERTION" && length($5)<50{LENGTH=length($5); print genome OFS "INSERTION" OFS LENGTH}$8=="DELETION" && length($4)<50{LENGTH=length($4); print genome OFS "DELETION" OFS LENGTH}' >> /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Comparative/Kestrel_Alignments/SharedVariation/Uniques/VariantCounts/UniqueKestrelAlignIndelSize.txt
done

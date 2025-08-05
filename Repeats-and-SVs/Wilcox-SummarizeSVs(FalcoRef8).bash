#Retrieving Summary of Unique Structural Variants
srun --pty -n 8 /bin/bash
rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/Summary
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/Summary
echo -e "Sample\tDEL_COUNT\tINS_COUNT\tINV_COUNT\tSegDEL_COUNT\tSegDUP_COUNT\tDEL_BP\tINS_BP\tINV_BP\tSegDEL_BP\tSegDUP_BP" > /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/Summary/strVariantSummary.txt

function summarize_SV {
Sample=$s
DEL_COUNT=`wc /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV.kestCoord.deletion50.uniq.bed | awk '{print $1}'` >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/Summary/strVariantSummary.txt
INS_COUNT=`wc /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV."$s"coord.insertion50.uniq.bed | awk '{print $1}'` >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/Summary/strVariantSummary.txt
INV_COUNT=`wc /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV."$s"coord.inversion.uniq.bed | awk '{print $1}'` >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/Summary/strVariantSummary.txt
SegDEL_COUNT=`wc /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV.kestCoord.seg_deletion.uniq.bed | awk '{print $1}'` >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/Summary/strVariantSummary.txt
SegDUP_COUNT=`wc /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV."$s"coord.seg_duplication.uniq.bed | awk '{print $1}'` >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/Summary/strVariantSummary.txt

DEL_BP=`awk 'BEGIN{basepairs=0}{size=$3-$2}{basepairs=basepairs+size}END{print basepairs}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV.kestCoord.deletion50.uniq.bed`
INS_BP=`awk 'BEGIN{basepairs=0}{size=$3-$2}{basepairs=basepairs+size}END{print basepairs}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV."$s"coord.insertion50.uniq.bed`
INV_BP=`awk 'BEGIN{basepairs=0}{size=$3-$2}{basepairs=basepairs+size}END{print basepairs}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV."$s"coord.inversion.uniq.bed`
SegDEL_BP=`awk 'BEGIN{basepairs=0}{size=$3-$2}{basepairs=basepairs+size}END{print basepairs}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV.kestCoord.seg_deletion.uniq.bed`
SegDUP_BP=`awk 'BEGIN{basepairs=0}{size=$3-$2}{basepairs=basepairs+size}END{print basepairs}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/SV/KestrelAlign/"$s"/Variants/Uniques/"$s"-v-Kestrel100KBhap.SV."$s"coord.seg_duplication.uniq.bed`

echo -e "$Sample\t$DEL_COUNT\t$INS_COUNT\t$INV_COUNT\t$SegDEL_COUNT\t$SegDUP_COUNT\t$DEL_BP\t$INS_BP\t$INV_BP\t$SegDEL_BP\t$SegDUP_BP" >> /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeat-SV/Summary/strVariantSummary.txt
}

for s in A D B FpBSx77 C E F Kestrel
do
summarize_SV &
done
wait

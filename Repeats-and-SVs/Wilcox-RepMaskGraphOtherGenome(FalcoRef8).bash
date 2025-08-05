#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32GB
#SBATCH --partition=cs
#SBATCH -o job.%J.out
#SBATCH --time=17:00:00

module purge
module load perl/intel/5.32.0

rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/DivCurves
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/DivCurves

function RepMaskGraph {
rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/"$o"/Analysis
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/"$o"/Analysis
cd /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/"$o"/Analysis
GENOME=`echo $o | perl -pe 's#Thrush#GCA_009819885.1_bCatUst1.pri_genomic.fna#' | perl -pe 's#Kakapo#GCF_004027225.2_bStrHab1.2.pri_genomic.fna#' | perl -pe 's#Seriema#GCA_009819825.1_bCarCri1.pri_genomic.fna#'`
GenomeSize=`grep -P "^total length:" /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/"$o"/"$GENOME".tbl | perl -pe "s#^total length:\s+([0-9]+) bp.+\([0-9]+.+#\1#"` 
Title=`echo "Interspersed Repeat Landscape: $o"`
perl /home/jw5478/Programs/RepeatMasker/util/calcDivergenceFromAlign.pl  -s /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/"$o"/Analysis/"$GENOME".align.divsum -a /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/"$o"/Analysis/"$GENOME".new.align /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/"$o"/"$GENOME".align
sed -E -i "s#;size=[0-9]+##g" /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/"$o"/Analysis/"$GENOME".align.divsum
perl /home/jw5478/Programs/RepeatMasker/util/createRepeatLandscape.pl -div /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/"$o"/Analysis/"$GENOME".align.divsum -t "$Title" -g $GenomeSize > /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/"$o"/Analysis/"$o".align.divsum.html
cp /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/"$o"/Analysis/"$o".align.divsum.html /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/OtherGenomes/rmBLAST/DivCurves/"$o".align.divsum.html
}

for o in Thrush Kakapo Seriema
do
RepMaskGraph &
done
wait

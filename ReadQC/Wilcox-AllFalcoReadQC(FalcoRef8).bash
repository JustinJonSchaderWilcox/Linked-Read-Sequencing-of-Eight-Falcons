#!/bin/csh
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=164:00:00
#SBATCH --mem=112G
# Output and error files
#SBATCH -o job.%J.out



module load gencore/1
module load jdk
module load jellyfish
module load gencore_qc/1.0

rm -r /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/
mkdir /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/
mkdir /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/FastQC
mkdir /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/FastQC/Merged
mkdir /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/Jellyfish

foreach s (A B C D E F FpBS Kestrel_Ft1)

echo "gunzip "$s" and remove 10X barcodes and merge raw reads"
rm -r /scratch/jw5478/FalconProject/RawData/RefGenomes/"$s"/MergedReads
mkdir /scratch/jw5478/FalconProject/RawData/RefGenomes/"$s"/MergedReads
gunzip -c /scratch/jw5478/FalconProject/RawData/RefGenomes/"$s"/"$s"_*R1*.fastq.gz | perl -pe 's#^[ATGCYRWSKMDVHBN]{26}##'  > /scratch/jw5478/FalconProject/RawData/RefGenomes/"$s"/MergedReads/"$s".mergedRR.fastq
gunzip -c /scratch/jw5478/FalconProject/RawData/RefGenomes/"$s"/"$s"_*R2*.fastq.gz >> /scratch/jw5478/FalconProject/RawData/RefGenomes/"$s"/MergedReads/"$s".mergedRR.fastq


echo "run FastQC on merged reads from "$s""


/project/jw5478/Programs/FastQC/fastqc /scratch/jw5478/FalconProject/RawData/RefGenomes/"$s"/MergedReads/"$s".mergedRR.fastq
mv /scratch/jw5478/FalconProject/RawData/RefGenomes/"$s"/MergedReads/*.html /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/FastQC/Merged
rm /scratch/jw5478/FalconProject/RawData/RefGenomes/"$s"/MergedReads/*.zip

echo "run Jellyfish on merged reads from "$s""

cat /scratch/jw5478/FalconProject/RawData/RefGenomes/"$s"/MergedReads/"$s".mergedRR.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}'  > /scratch/jw5478/FalconProject/RawData/RefGenomes/"$s"/MergedReads/"$s".mergedRR.fasta

rm -r /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/Jellyfish/"$s"
mkdir /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/Jellyfish/"$s"
cd /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/Jellyfish/"$s"
jellyfish count -m 27 -s 100M -t 28 -o /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/Jellyfish/"$s"/"$s".mergedRR.jf -C /scratch/jw5478/FalconProject/RawData/RefGenomes/"$s"/MergedReads/"$s".mergedRR.fasta
jellyfish histo -o "$s".mergedRR.jf.histo /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/Jellyfish/"$s"/"$s".mergedRR.jf
jellyfish stats -v --output="/scratch/jw5478/FalconProject/RawData/RefGenomes/QC/Jellyfish/"$s"/"$s"_stats.output.jf" /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/Jellyfish/"$s"/"$s".mergedRR.jf


echo "run kat on merged reads from "$s""
setenv LD_LIBRARY_PATH /scratch/gencore/.local/easybuild/software/gencore_anaconda2/4.2.0/pkgs/boost-1.61.0-py27_1/lib/
/scratch/gencore/.local/easybuild/software/gencore_anaconda2/4.2.0/pkgs/kat-2.3.1-boost1.61_2/bin/kat gcp -t 28 -o "$s"_mergedRR /scratch/jw5478/FalconProject/RawData/RefGenomes/"$s"/MergedReads/"$s".mergedRR.fasta

set Title=`grep -P "# Title:" /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/Jellyfish/"$s"/"$s"_mergedRR.mx | perl -pe "s\#\\" | sed "s\ Title:\\"`
set Xlab=`grep -P "# XLabel:" /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/Jellyfish/"$s"/"$s"_mergedRR.mx | perl -pe "s\#\\" | sed "s\ XLabel:\\"`
set Ylab=`grep -P "# YLabel:" /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/Jellyfish/"$s"/"$s"_mergedRR.mx | perl -pe "s\#\\" | sed "s\ YLabel:\\"`
/scratch/gencore/.local/easybuild/software/gencore_anaconda2/4.2.0/pkgs/kat-2.3.1-boost1.61_2/bin/kat plot density -p png -t "$Title" -a "$Xlab" -b "$Ylab" -o "$s"_mergedRR.GC_density.png /scratch/jw5478/FalconProject/RawData/RefGenomes/QC/Jellyfish/"$s"/"$s"_mergedRR.mx

end

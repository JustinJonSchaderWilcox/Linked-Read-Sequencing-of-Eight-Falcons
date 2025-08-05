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
module load gencore/1
module load gencore_qc/1.0
module load gencore_variant_detection/1.0
module load jdk

#Variables
s=$1
FalcoMito="/scratch/jw5478/FalconProject/Projects/RefGenome8/Mitochondria/Mitogenome"
longranger="/scratch/gencore/software/longranger/longranger-2.2.2/longranger"
TYPE=`echo "$s" | perl -pe "s#\bB\b#Peregrine#" | perl -pe "s#\bC\b#Peregrine#" | perl -pe "s#\bFpBSx77\b#Peregrine#" | perl -pe "s#\bA\b#Hierofalco#" | perl -pe "s#\bD\b#Hierofalco#" | perl -pe "s#\bE\b#Hierofalco#" | perl -pe "s#\bF\b#Hierofalco#"`
DATE=`ls $FalcoMito/Downloads/Sequences/Alignments/Consensus/"$TYPE"MitoGenome_*.cons.fasta | perl -pe "s#.*"$TYPE"MitoGenome_(.*)\.cons.fasta#\1#"`



cd $FalcoMito/LongRanger/WGS
r=`echo $s | perl -pe 's#FpBSx77#FpBS#' | perl -pe 's#Kestrel#Kestrel_Ft1#'`
$longranger wgs \
--id="$s"_MitoLR \
--fastqs=/scratch/jw5478/RawData/RefGenomes/"$r"/ \
--indices=SI-GA-D1 \
--reference=$FalcoMito/LongRanger/ConMitoRef/refdata-"$TYPE"MitoGenome_"$DATE".cons/ \
--vcmode=gatk:/scratch/gencore/.local/easybuild/software/gencore_variant_detection/1.0/opt/gatk-3.5/GenomeAnalysisTK.jar \
--sample="$r" \
--sex="f" \
--localcores=28

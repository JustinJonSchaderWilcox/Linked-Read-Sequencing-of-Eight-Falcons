#!/bin/csh
#SBATCH --mem=116GB
#SBATCH -p serial
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=28
# Output and error files
#SBATCH -o job.%J.out

#This program takes an input sample name and then produces a phased scaffolds in fasta format, while excluding all phaseblocks less than 16KB


set supernova="/scratch/gencore/supernova/supernova-2.1.1/supernova"

cd /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2

$supernova mkoutput \
--style=pseudohap2 \
--asmdir=/scratch/gencore/novaseq/190513_A00534_0032_AHF3WFDSXX/FpBSx77/outs/assembly \
--outprefix=FpBSx77_16KBph2 \
--minsize=16000 \
--headers=short


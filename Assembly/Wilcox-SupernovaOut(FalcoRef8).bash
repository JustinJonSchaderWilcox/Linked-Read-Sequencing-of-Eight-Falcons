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

#This program takes an input sample name and then produces three different mkoutput from 10X supernova
set s=$1

set supernova="/scratch/gencore/supernova/supernova-2.1.1/supernova"

cd /scratch/jw5478/FalconProject/Projects/RefGenome6/SuperNovaOut

$supernova mkoutput \
--style=raw \
--asmdir=/scratch/gencore/novaseq/181220_A00534_0025_AHHM57DMXX/"$s"/outs/assembly \
--outprefix="$s"_raw2 \
--headers=short

$supernova mkoutput \
--style=megabubbles \
--asmdir=/scratch/gencore/novaseq/181220_A00534_0025_AHHM57DMXX/"$s"/outs/assembly \
--outprefix="$s"_16KBmegabubbles \
--minsize=16000 \
--headers=short

$supernova mkoutput \
--style=pseudohap2 \
--asmdir=/scratch/gencore/novaseq/181220_A00534_0025_AHHM57DMXX/"$s"/outs/assembly \
--outprefix="$s"_16KBph2 \
--minsize=16000 \
--headers=short


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

s=$1
/project/jw5478/Programs/psmc/utils/splitfa /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/ConsensusGenomes/"$s"_16KBph2.haploidConsensus.psmcfa > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/ConsensusGenomes/"$s"_16KBph2.haploidConsensus.split.psmcfa

for b in $(seq 1 100)
do
/project/jw5478/Programs/psmc/psmc -N100 -t 12 -r 2 -b -p "8+28*3+8" -o /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/BootStrap/Rounds/round-"$b"."$s"_16KBph2.diploid.psmc /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/ConsensusGenomes/"$s"_16KBph2.haploidConsensus.split.psmcfa &
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 28 ]
do
sleep 1m
JOBS=`jobs | wc | awk '{print $1}'`
done
done
wait
cat /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/"$s"_16KBph2.diploid.psmc /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/BootStrap/Rounds/round*"$s"_16KBph2.diploid.psmc > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/BootStrap/"$s"_16KBph2.diploid.combined.psmc

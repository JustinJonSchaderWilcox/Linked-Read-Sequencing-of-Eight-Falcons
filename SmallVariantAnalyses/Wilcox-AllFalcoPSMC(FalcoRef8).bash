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
module load gencore/1
module load gencore_variant_detection/1.0
module load perl

samtools mpileup -uf /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta /scratch/jw5478/FalconProject/Projects/RefGenome8/LongrangerOut/SelfAlign/16KBph2Haploid/"$s"_16KBph2-HaploidLR/outs/phased_possorted_bam.bam | bcftools call -c --threads 28  - | bcftools filter --SnpGap 10 - | bcftools view  -i  'MIN(DP)>20 & MAX(DP)<120'   - > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/ConsensusGenomes/Temp/"$s"_16KBph2.100KB.haploid.allsite.vcf

vcfutils.pl vcf2fq /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/ConsensusGenomes/Temp/"$s"_16KBph2.100KB.haploid.allsite.vcf > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/ConsensusGenomes/"$s"_diploid.fq
gzip /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/ConsensusGenomes/"$s"_diploid.fq
/project/jw5478/Programs/psmc/utils/fq2psmcfa -q 30 /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/ConsensusGenomes/"$s"_diploid.fq.gz > /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/ConsensusGenomes/"$s"_16KBph2.haploidConsensus.psmcfa
/project/jw5478/Programs/psmc/psmc -N100 -t 12 -r 2 -p "8+28*3+8" -o /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/"$s"_16KBph2.diploid.psmc /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Demography/PSMC/ConsensusGenomes/"$s"_16KBph2.haploidConsensus.psmcfa

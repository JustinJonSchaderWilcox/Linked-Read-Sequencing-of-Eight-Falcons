#!/bin/csh
#SBATCH --partition=condo
#SBATCH --mem=460GB
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=28
# Output and error files
#SBATCH -o job.%J.out

module purge
module load gencore/1
module load gencore_qc/1.0
module load gencore_variant_detection/1.0
module load jdk


set longranger="/scratch/gencore/software/longranger/longranger-2.2.2/longranger"

set s=$1
set SEX=$2



rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques
grep -P ">[0-9]+_[12]" /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_16KBph2.diploid.fasta | perl -pe "s#(>[0-9]+_[12]).*#\1#" > /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_AllScaf.text
grep -P ">[0-9]+_[1]" /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_16KBph2.diploid.fasta | perl -pe "s#(>[0-9]+_[1]).*#\1#" > /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS1Scaf.text
grep -P ">[0-9]+_[2]" /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_16KBph2.diploid.fasta | perl -pe "s#(>[0-9]+_[2]).*#\1#" > /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS2Scaf.text

cat /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS1Scaf.text > /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS1Scaf.uniq-int.text
cat /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS2Scaf.text > /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS2Scaf.uniq-int.text

set PH2_Names=`cat /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS2Scaf.text`
foreach n ($PH2_Names)
set r=`echo $n | sed 's#_2#_1#'`
sed -i "s#$r##" /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS1Scaf.uniq-int.text
end
perl -pe -chomp /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS1Scaf.uniq-int.text | perl -pe 's#_1#_1\n#g' > /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS1Scaf.uniq.text
rm /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS1Scaf.uniq-int.text

set PH1_Names=`cat /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS1Scaf.text`
foreach n ($PH1_Names)
set r=`echo $n | sed 's#_1#_2#'`
sed -i "s#$r##" /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS2Scaf.uniq-int.text
end
perl -pe -chomp /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS2Scaf.uniq-int.text | perl -pe 's#_2#_2\n#g' > /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS2Scaf.uniq.text
rm /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS2Scaf.uniq-int.text

cat /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS1Scaf.text /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_PS2Scaf.uniq.text > /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_IncHapRefScaf.text


cat /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_16KBph2.diploid.fasta | perl -pe "s#^([A-Z]{80})\n#\1#" | perl -pe "s#([A-Z]{80,})(>.+\n)#\1\n\2#" > /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/DiploidGenomesMod/"$s"_16KBph2.diploid.sl.fasta

set hapScaffolds=`cat /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/"$s"_LG-PS_Uniques/"$s"_IncHapRefScaf.text`

rm /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/Haploid16KBph2Seed/"$s"_16KBph2.haploid.fasta
foreach h ($hapScaffolds)
grep -A 1 "$h" /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/DiploidGenomesMod/"$s"_16KBph2.diploid.sl.fasta >> /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/Haploid16KBph2Seed/"$s"_16KBph2.haploid.fasta
end

cd /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/
rm /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/DiploidGenomes/100kbDiploid/"$s"_16KBph2.100KB.haploid.fasta
fastaq filter --min_length 100000 /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/Haploid16KBph2Seed/"$s"_16KBph2.haploid.fasta /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta


$longranger mkref /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/"$s"_16KBph2.100KB.haploid.fasta

picard CreateSequenceDictionary R=/scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/refdata-"$s"_16KBph2.100KB.haploid/fasta/genome.fa O=/scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/refdata-"$s"_16KBph2.100KB.haploid/fasta/genome.dict

cd /scratch/jw5478/FalconProject/Projects/RefGenome8/LongrangerOut/SelfAlign/16KBph2Haploid

$longranger wgs \
--id="$s"_16KBph2-HaploidLR \
--fastqs=/scratch/jw5478/FalconProject/RawData/RefGenomes/Kestrel_Ft1/ \
--indices=SI-GA-D1 \
--reference=/scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/refdata-"$s"_16KBph2.100KB.haploid \
--vcmode=gatk:/scratch/gencore/.local/easybuild/software/gencore_variant_detection/1.0/opt/gatk-3.5/GenomeAnalysisTK.jar \
--sample=Kestrel_Ft1 \
--sex="$SEX" \
--localcores=28


rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/LongrangerOut/SelfAlign/16KBph2Haploid/"$s"_16KBph2-HaploidLR/PHASER_SVCALLER_CS/
rm /scratch/jw5478/FalconProject/Projects/RefGenome8/LongrangerOut/SelfAlign/16KBph2Haploid/"$s"_16KBph2-HaploidLR/_*



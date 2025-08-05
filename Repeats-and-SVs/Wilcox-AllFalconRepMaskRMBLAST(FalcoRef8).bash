#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180GB
#SBATCH -o job.%J.out
#SBATCH --time=168:00:00

#Modules
module purge
module load perl/intel/5.32.0

#Input
s=$1

##RepeatMasker
rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/rmBLAST/"$s"
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/rmBLAST/"$s"

cd /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/rmBLAST/"$s"
/home/jw5478/Programs/RepeatMasker/RepeatMasker -e rmblast -a -lib /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/All_FalcoRefGenomeRepModDB-ID80families.annotated.fa -dir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/rmBLAST/"$s" -gccalc -poly -s -pa 12 /scratch/jw5478/FalconProject/Projects/RefGenomes/Genomes/SuperNova10X/16KBph2Diploid/"$s"_16KBph2.diploid.fasta

rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/rmBLAST/"$s"/Analysis
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/rmBLAST/"$s"/Analysis
cd /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/rmBLAST/"$s"/Analysis

GenomeSize=`grep -P "^total length:" /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/rmBLAST/"$s"/"$s"_16KBph2.diploid.fasta.tbl | perl -pe "s#^total length:\s+([0-9]+) bp.+\([0-9]+.+#\1#"` 
Title=`grep -P "^$s\t" /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/SampleKey.txt | perl -pe "s#^[a-zA-Z0-9]+\s+([A-Z][a-zA-Z()\s]+).*#Interspersed Repeat Landscape: \1#"`
perl /home/jw5478/Programs/RepeatMasker/util/calcDivergenceFromAlign.pl  -s /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/rmBLAST/"$s"/Analysis/"$s"_16KBph2.diploid.fasta.align.divsum -a /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/rmBLAST/"$s"/Analysis/"$s"_16KBph2.diploid.fasta.new.align /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/rmBLAST/"$s"/"$s"_16KBph2.diploid.fasta.align
perl /home/jw5478/Programs/RepeatMasker/util/createRepeatLandscape.pl -div /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/rmBLAST/"$s"/Analysis/"$s"_16KBph2.diploid.fasta.align.divsum -t "$Title" -g $GenomeSize > /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/rmBLAST/"$s"/Analysis/"$s"_16KBph2.diploid.fasta.align.divsum.html

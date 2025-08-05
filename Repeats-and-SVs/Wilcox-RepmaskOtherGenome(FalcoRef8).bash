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
rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"

rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/rmBLAST
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/rmBLAST

cd /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/rmBLAST
/home/jw5478/Programs/RepeatMasker/RepeatMasker -e rmblast -a -lib /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/All_FalcoRefGenomeRepModDB-ID80families.annotated.fa -dir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/rmBLAST -gccalc -poly -s -pa 12 /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/Genomes/"$s".fasta

rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/rmBLAST/Analysis
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/rmBLAST/Analysis
cd /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/rmBLAST/Analysis

GenomeSize=`grep -P "^total length:" /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/rmBLAST/"$s".fasta.tbl | perl -pe "s#^total length:\s+([0-9]+) bp.+\([0-9]+.+#\1#"` 
Title=`grep -P "^$s\t" /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomeKey.txt | perl -pe "s#^[a-zA-Z0-9_]+\s+([A-Z][a-zA-Z()\s]+).*#Interspersed Repeat Landscape: \1#"`
perl /home/jw5478/Programs/RepeatMasker/util/calcDivergenceFromAlign.pl  -s /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/rmBLAST/Analysis/"$s".fasta.align.divsum -a /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/rmBLAST/Analysis/"$s".fasta.new.align /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/rmBLAST/"$s".fasta.align
perl /home/jw5478/Programs/RepeatMasker/util/createRepeatLandscape.pl -div /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/rmBLAST/Analysis/"$s".fasta.align.divsum -t "$Title" -g $GenomeSize > /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/rmBLAST/Analysis/"$s".fasta.align.divsum.html

##RepeatMasker
rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/HMMER
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/HMMER

cd /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/HMMER
/home/jw5478/Programs/RepeatMasker/RepeatMasker -e hmmer -a -lib /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/All_FalcoRefGenomeRepModDB-families.hmm -dir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/HMMER -gccalc -poly -s -pa 24 /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/Genomes/"$s".fasta

rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/HMMER/Analysis
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/HMMER/Analysis
cd /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/HMMER/Analysis

GenomeSize=`grep -P "^total length:" /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/HMMER/"$s".fasta.tbl | perl -pe "s#^total length:\s+([0-9]+) bp.+\([0-9]+.+#\1#"` 
Title=`grep -P "^$s\t" /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomeKey.txt | perl -pe "s#^[a-zA-Z0-9_]+\s+([A-Z][a-zA-Z()\s]+).*#Interspersed Repeat Landscape: \1#"`
perl /home/jw5478/Programs/RepeatMasker/util/calcDivergenceFromAlign.pl  -s /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/HMMER/Analysis/"$s".fasta.align.divsum -a /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/HMMER/Analysis/"$s".fasta.new.align /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/HMMER/"$s".fasta.align
perl /home/jw5478/Programs/RepeatMasker/util/createRepeatLandscape.pl -div /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/HMMER/Analysis/"$s".fasta.align.divsum -t "$Title" -g $GenomeSize > /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/RepeatMasker/OtherGenomes/"$s"/HMMER/Analysis/"$s".fasta.align.divsum.html

#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --mem=118GB
#SBATCH -p serial
# Set number of nodes to run
# Set number of tasks to run
#SBATCH --ntasks=28
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

module purge
module load perl
module load gcc

for r in LINE LTR DNA
do
/home/jw5478/Executables/Programs/Mafft/mafft-linsi --nuc --reorder --thread 28 --quiet --maxiterate 1000 /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Gene_Search/Sequences/"$r"_Sequences.fasta > /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Gene_Search/Alignments/"$r"_Sequences.linsi.fasta
done

for r in LINE LTR DNA
do
/home/jw5478/Executables/Programs/HMMER/bin/hmmbuild --cpu 28 /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Gene_Search/Alignments/"$r"_Sequences.linsi.hmm /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Gene_Search/Alignments/"$r"_Sequences.linsi.fasta
done

for s in A B C D E F FpBSx77 Kestrel
do
for r in LINE LTR DNA
do
TARGET_SIZE=`echo $r | perl -pe 's#LINE#2000#' | perl -pe 's#LTR#2200#' | perl -pe 's#DNA#1200#'`
/home/jw5478/Executables/Programs/HMMER/bin/nhmmer --cpu 28 -A /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Gene_Search/ExtractedSequences/"$s"_extracted"$r".stk -E 0.000001 --incE 0.000001 /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Gene_Search/Alignments/"$r"_Sequences.linsi.hmm /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Active/ORFs/"$s"_16KBph2.diploid.fasta.all_rep.ORF.fasta
/home/jw5478/Executables/Programs/HMMER/bin/esl-reformat fasta /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Gene_Search/ExtractedSequences/"$s"_extracted"$r".stk | awk 'NR==1{print $0; next}/^>/{print "\n"$0;next}!/^>/{ printf($0)}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Gene_Search/ExtractedSequences/"$s"_extracted"$r".fasta
awk -v target_size=$TARGET_SIZE '/>/{TE_Name=$0; next}!/>/{ORF_SIZE=length($0)}ORF_SIZE>=target_size{print TE_Name";""length""="ORF_SIZE; print $0}' /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Gene_Search/ExtractedSequences/"$s"_extracted"$r".fasta > /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Gene_Search/CandidateActiveORF/"$s"_candidate"$r".fasta
/home/jw5478/Executables/Programs/USEARCH/usearch11.0.667_i86linux32 -search_local /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Gene_Search/ExtractedSequences/"$s"_extracted"$r".fasta -id -0.75 -strand both -target_cov 0.9 -threads 28 -db /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Gene_Search/Sequences/"$r"_Sequences.fasta -matched  /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepeatMasker/Ultimate/Gene_Search/AnnotatedActiveORF/"$s"_annotated"$r".fasta
done
done

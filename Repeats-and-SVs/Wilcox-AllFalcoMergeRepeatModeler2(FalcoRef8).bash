#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180GB
#SBATCH -o job.%J.out
#SBATCH --time=168:00:00

#Script to compile Repeats from all Stockholm files into a single file

rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp
for s in A B C D E F Kestrel FpBSx77
do 
perl -pe "s#>#>$s\_#" /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/Single_Falco/"$s"_FalcoRefGenomeRepModDB-families.fa | perl -pe 's#\s\(.*Final Multiple Alignment (Size = [0-9]+).*#;$1#' | perl -pe 's#Size = #size=#' >> /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-families.fa
done

/home/jw5478/Programs/USEARCH/usearch11.0.667_i86linux32 -sortbysize /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-families.fa -fastaout /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-families.size_filt.fa -minsize 20
/home/jw5478/Programs/USEARCH/usearch11.0.667_i86linux32 -sortbylength /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-families.size_filt.fa -fastaout /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-families.size_filter.sort_length.fa 
/home/jw5478/Programs/USEARCH/usearch11.0.667_i86linux32 -cluster_fast /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-families.size_filter.sort_length.fa -threads 48 -id 0.8 -sort length -maxaccepts 1 -maxrejects 0 -strand both -sizein -sizeout -consout /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-ID80families.fa -uc /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-ID80families.uc

/home/jw5478/Programs/USEARCH/usearch11.0.667_i86linux32 -sortbysize /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-ID80families.fa -fastaout /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-ID80families.size_filter.fa -minsize 80

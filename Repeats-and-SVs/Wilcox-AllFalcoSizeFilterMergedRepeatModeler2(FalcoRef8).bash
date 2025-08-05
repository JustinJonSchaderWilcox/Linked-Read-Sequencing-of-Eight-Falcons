srun --pty -c 4 /bin/bash
rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco
mkdir /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco
cp -r /home/jw5478/Temp /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco
/home/jw5478/Programs/USEARCH/usearch11.0.667_i86linux32 -sortbysize /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-ID80families.fa -fastaout /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-ID80families.size_filter40.fa -minsize 40

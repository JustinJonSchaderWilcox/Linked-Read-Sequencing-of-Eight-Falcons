#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180GB
#SBATCH -o job.%J.out
#SBATCH --time=168:00:00

module purge
rm -r /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp
cp -r /home/jw5478/Temp /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco
rm /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/All_FalcoRefGenomeRepModDB-ID80families.annotated.fa
/home/jw5478/Programs/USEARCH/usearch11.0.667_i86linux32 -sortbysize /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-ID80families.fa -fastaout /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-ID80families.size_filter40.fa -minsize 40
cat /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-ID80families.size_filter40.fa > /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/All_FalcoRefGenomeRepModDB-ID80families.annotated.fa
for f in $(grep ">" /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-ID80families.size_filter40.fa | perl -pe 's#>##')
do
CLUSTER_NUM=`echo $f | perl -pe 's#;.*##' | perl -pe 's#Cluster##'`
i=`awk -v cluster=$CLUSTER_NUM '$2==cluster && $10=="*"{print $9; exit}' /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/Temp/All_FalcoRefGenomeRepModDB-ID80families.uc | perl -pe 's#;.*##'`
sed -i "s%$f\$%$i%" /scratch/jw5478/FalconProject/Projects/RefGenomes/Repeats/16KBph2Diploid/All_Falco/All_FalcoRefGenomeRepModDB-ID80families.annotated.fa
done
rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/Repeats/16KBph2Diploid/RepModDB/All_Falco/Temp

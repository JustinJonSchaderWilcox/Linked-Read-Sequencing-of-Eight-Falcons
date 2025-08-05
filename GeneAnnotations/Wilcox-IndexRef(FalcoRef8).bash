#Indexing Reference Genome
module purge
module load gencore/1
module load gencore_annotation/1.0

for s in A D B FpBSx77 C E F Kestrel
do
rm /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/Haploid16KBph2Seed/"$s"_16KBph2.haploid.fasta.fai
samtools faidx /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/Haploid16KBph2Seed/"$s"_16KBph2.haploid.fasta
done

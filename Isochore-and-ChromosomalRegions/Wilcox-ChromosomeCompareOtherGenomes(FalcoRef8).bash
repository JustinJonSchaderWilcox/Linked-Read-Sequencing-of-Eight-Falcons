#Downloading genomes of non-falconidae Inopinaves for annotatioon of past chromosome states 
#Downloading Chromosomal Scale Scaffolds for Swainson's Thrush
##chose Swainson's Thrush genome as it was the only genome that was chromosome assembled and had a normal avian karyotype; based on search of txid9126[Organism:exp] and txid9224[Organism:exp] for chromosome assemblies on NCBI Genbank https://www.ncbi.nlm.nih.gov/
##2N equal 80 plus sex chromosomes
rm -r /scratch/jw5478/FalconProject/OtherGenomes/
mkdir /scratch/jw5478/FalconProject/OtherGenomes/
rm -r /scratch/jw5478/FalconProject/OtherGenomes/Catharus_ustulatus
mkdir /scratch/jw5478/FalconProject/OtherGenomes/Catharus_ustulatus
cd /scratch/jw5478/FalconProject/OtherGenomes/Catharus_ustulatus
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Catharus_ustulatus/latest_assembly_versions/GCF_009819885.1_bCatUst1.pri/GCF_009819885.1_bCatUst1.pri_genomic.fna.gz
gunzip /scratch/jw5478/FalconProject/OtherGenomes/Catharus_ustulatus/GCF_009819885.1_bCatUst1.pri_genomic.fna.gz
#Downloading Chromosomal Scale Scaffolds for Kakapo
##chose Kakapo (Strigops habroptila) as it has as is only chromosome scale assembly from a  parrot, and has low chromosome count; txid9224[Organism:exp]  on https://www.ncbi.nlm.nih.gov/
##2N=23 plus sex chromosomes
rm -r /scratch/jw5478/FalconProject/OtherGenomes/Strigops_habroptila
mkdir /scratch/jw5478/FalconProject/OtherGenomes/Strigops_habroptila
cd /scratch/jw5478/FalconProject/OtherGenomes/Strigops_habroptila
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/027/225/GCF_004027225.2_bStrHab1.2.pri/GCF_004027225.2_bStrHab1.2.pri_genomic.fna.gz
gunzip /scratch/jw5478/FalconProject/OtherGenomes/Strigops_habroptila/GCF_004027225.2_bStrHab1.2.pri_genomic.fna.gz
#Downloading Chromosomal Scale Scaffolds for Red-legged seriema (Cariama cristata)
##chose Cariama cristata as it has as is only chromosome scale assembly for Cariamiformes  (seriamas), and has high chromosome count; txid1956187[Organism:exp] on https://www.ncbi.nlm.nih.gov/
##2N=100 plus sex chromosomes
rm -r /scratch/jw5478/FalconProject/OtherGenomes/Cariama_cristata
mkdir /scratch/jw5478/FalconProject/OtherGenomes/Cariama_cristata
cd /scratch/jw5478/FalconProject/OtherGenomes/Cariama_cristata
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/819/825/GCA_009819825.1_bCarCri1.pri/GCA_009819825.1_bCarCri1.pri_genomic.fna.gz
gunzip /scratch/jw5478/FalconProject/OtherGenomes/Cariama_cristata/GCA_009819825.1_bCarCri1.pri_genomic.fna.gz

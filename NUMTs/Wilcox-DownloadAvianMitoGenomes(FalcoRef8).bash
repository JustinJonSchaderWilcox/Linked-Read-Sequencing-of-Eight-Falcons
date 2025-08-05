#Script to download all Falconidae Mitogenomes with E-Utilities
#Download All Falcondiae Mitogenomes with E-Utilities
module purge
module load perl
FalconidaeMito="/scratch/jw5478/FalconProject/Projects/RefGenome8/Mitochondria/NUMT"
edirect="/home/jw5478/Executables/Programs/edirect"
DATE=`date '+ %d%B%y' | perl -pe 's#\s##g'`
#Download Falcondiae Mitogenomes
$edirect/esearch -db nucleotide -query "txid8949[Organism] mitochondrion[filter] AND ("16000"[SLEN] : "20000"[SLEN])" > $FalconidaeMito/Downloads/Records/FalconidaeMitoGenome_"$DATE"_DownloadRecord.txt
$edirect/esearch -db nucleotide -query "txid8949[Organism] mitochondrion[filter] AND ("16000"[SLEN] : "20000"[SLEN])" | $edirect/efetch -db nucleotide -format fasta >  $FalconidaeMito/Downloads/Sequences/FalconidaeMitoGenome_"$DATE".fasta
#Download other Mitogenomes
#Swainson's Thrush
$edirect/esearch -db nucleotide -query "CM020378 AND txid91951[Organism] mitochondrion[filter] AND ("16000"[SLEN] : "20000"[SLEN])" | $edirect/efetch -db nucleotide -format fasta >  $FalconidaeMito/Downloads/Sequences/ThrushNonFalconidaeMitoGenomeCM020378.fasta
#Cacatua alba
$edirect/esearch -db nucleotide -query "MT920475 AND txid9223[Organism] mitochondrion[filter] AND ("16000"[SLEN] : "20000"[SLEN])" | $edirect/efetch -db nucleotide -format fasta >  $FalconidaeMito/Downloads/Sequences/CockatooNonFalconidaeMitoGenomeAP010797.fasta
#Picus canis
$edirect/esearch -db nucleotide -query "NC_045372.1 AND txid9220[Organism] mitochondrion[filter] AND ("16000"[SLEN] : "20000"[SLEN])" | $edirect/efetch -db nucleotide -format fasta >  $FalconidaeMito/Downloads/Sequences/PicusNonFalconidaeMitoGenomeAP010797.fasta
#Strix occidentalis
$edirect/esearch -db nucleotide -query "MF431746 AND txid30458[Organism] mitochondrion[filter] AND ("16000"[SLEN] : "20000"[SLEN])" | $edirect/efetch -db nucleotide -format fasta >  $FalconidaeMito/Downloads/Sequences/OwlNonFalconidaeMitoGenomeAP010797.fasta
#Goshawk
$edirect/esearch -db nucleotide -query "AP010797 AND txid8957[Organism] mitochondrion[filter] AND ("16000"[SLEN] : "20000"[SLEN])" | $edirect/efetch -db nucleotide -format fasta >  $FalconidaeMito/Downloads/Sequences/GoshawkNonFalconidaeMitoGenomeAP010797.fasta
#Gavia stellata
$edirect/esearch -db nucleotide -query "AY293618 AND txid37038[Organism] mitochondrion[filter] AND ("16000"[SLEN] : "20000"[SLEN])" | $edirect/efetch -db nucleotide -format fasta >  $FalconidaeMito/Downloads/Sequences/LoonNonFalconidaeMitoGenomeAP010797.fasta
#Balearica regulorum
$edirect/esearch -db nucleotide -query "FJ769841 AND txid925459[Organism] mitochondrion[filter] AND ("16000"[SLEN] : "20000"[SLEN])" | $edirect/efetch -db nucleotide -format fasta >  $FalconidaeMito/Downloads/Sequences/CraneNonFalconidaeMitoGenomeAP010797.fasta
#Columba livia
$edirect/esearch -db nucleotide -query "KP319029 AND txid8932[Organism] mitochondrion[filter] AND ("16000"[SLEN] : "20000"[SLEN])" | $edirect/efetch -db nucleotide -format fasta >  $FalconidaeMito/Downloads/Sequences/DoveNonFalconidaeMitoGenomeAP010797.fasta
#Calypte anna (Humming bird)
$edirect/esearch -db nucleotide -query "CM016612 AND txid9242[Organism] mitochondrion[filter] AND ("16000"[SLEN] : "20000"[SLEN])" | $edirect/efetch -db nucleotide -format fasta >  $FalconidaeMito/Downloads/Sequences/HummingBirdNonFalconidaeMitoGenomeAP010797.fasta
#Gallus Gallus
$edirect/esearch -db nucleotide -query "AP003322 AND txid9031[Organism] mitochondrion[filter] AND ("16000"[SLEN] : "20000"[SLEN])" | $edirect/efetch -db nucleotide -format fasta >  $FalconidaeMito/Downloads/Sequences/GallusNonFalconidaeMitoGenomeAP003322.fasta
#Ostrich
$edirect/esearch -db nucleotide -query "AF338715 AND txid8801[Organism] mitochondrion[filter] AND ("16000"[SLEN] : "20000"[SLEN])" | $edirect/efetch -db nucleotide -format fasta >  $FalconidaeMito/Downloads/Sequences/StruthioNonFalconidaeMitoGenomeAF338715.fasta
cat $FalconidaeMito/Downloads/Sequences/*FalconidaeMitoGenome*.fasta > $FalconidaeMito/Downloads/Sequences/Compiled/AvesSampledFalconidaeMitoGenome_"$DATE".fasta

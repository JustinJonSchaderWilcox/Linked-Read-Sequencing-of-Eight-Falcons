#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --mem=118GB
#SBATCH -p serial
# Set number of tasks to run
#SBATCH --ntasks=28
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

module purge
module load perl
module load gencore/1
module load gencore_variant_detection/1.0
module load jdk
module load blast/2.3.0

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/Transcripts
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/ORF
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Merged
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/Modified
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Databases
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/BLAST_Outs

for s in A D B FpBSx77 C E F Kestrel
do
/scratch/gencore/.local/easybuild/software/gencore_annotation/1.0/opt/transdecoder/util/cufflinks_gtf_genome_to_cdna_fasta.pl /scratch/jw5478/FalconProject/Projects/RefGenome8/Transciptome-to-GeneAnnotation/Stringtie_GTF/"$s"_stringtie_merge.gtf /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/Haploid16KBph2Seed/"$s"_16KBph2.haploid.fasta > /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/Transcripts/"$s"_TranscriptSeqs.fasta &
/scratch/gencore/.local/easybuild/software/gencore_annotation/1.0/opt/transdecoder/util/cufflinks_gtf_to_alignment_gff3.pl /scratch/jw5478/FalconProject/Projects/RefGenome8/Transciptome-to-GeneAnnotation/Stringtie_GTF/"$s"_stringtie_merge.gtf > /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/Transcripts/"$s"_stringtie_merge.gff3 &
done
wait

cd /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/ORF
for s in A D B FpBSx77 C E F Kestrel
do
/scratch/gencore/.local/easybuild/software/gencore_annotation/1.0/opt/transdecoder/TransDecoder.LongOrfs -t /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/Transcripts/"$s"_TranscriptSeqs.fasta &
done
wait

for s in A D B FpBSx77 C E F Kestrel
do
/scratch/gencore/.local/easybuild/software/gencore_annotation/1.0/opt/transdecoder/util/cdna_alignment_orf_to_genome_orf.pl /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/ORF/"$s"_TranscriptSeqs.fasta.transdecoder_dir/longest_orfs.gff3 /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/Transcripts/"$s"_stringtie_merge.gff3 /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/Transcripts/"$s"_TranscriptSeqs.fasta > /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/ORF/"$s"_TranscriptSeqs.fasta.transdecoder.genome.gff3 &
done
wait
rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/Transcripts


cd /scratch/jw5478/FalconProject/Projects/RefGenome8/Scripts/EvolutionaryPressures/Submissions
#Run perl script to extract DNA & protein sequences for each gene for each species
for s in A D B FpBSx77 C E F Kestrel
do
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/"$s"
perl /scratch/ieh211/Ptychadena/dnds/paml/process_gff_cds_proteins.pl -g /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/ORF/"$s"_TranscriptSeqs.fasta.transdecoder.genome.gff3 -d /scratch/jw5478/FalconProject/Projects/RefGenome8/SuperNovaOut/16KBph2/16KBph2Mod/PH2_100kbRef/Haploid16KBph2Seed/"$s"_16KBph2.haploid.fasta -o /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/"$s" &
done
wait
rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/ORF

for s in A D B FpBSx77 C E F Kestrel
do
rm /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/"$s"/*.prot.fasta &
done
wait

#Rename gene-ids to include species name for all 
GENES=`ls /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/*/ | perl -pe 's#.*\/##' | perl -pe 's#.*-##' | grep ".dna." | sort -k1,1V | uniq`
for s in A D B FpBSx77 C E F Kestrel
do
for i in $GENES
do
N=`echo $i | perl -pe 's#.dna.fasta##' | perl -pe 's#MSTRG##'`
d=`expr $N / 1000`
mkdir -p /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/"$s"/"$d"
mv /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/"$s"/"$i" /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/"$s"/"$d"/"$s"-"$i" &
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 28 ]
do
sleep 1s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
done
wait

#Rename fasta headers to include species names
for s in A D B FpBSx77 C E F Kestrel
do
for i in $GENES
do
N=`echo $i | perl -pe 's#.dna.fasta##' | perl -pe 's#MSTRG##'`
d=`expr $N / 1000`
sed -i "s#>#>$s\_#" /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/"$s"/"$d"/"$s"-"$i" &
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 28 ]
do
sleep 1s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
done
wait

for s in A D B FpBSx77 C E F Kestrel
do
cat /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/*/*/$s-*.dna.fasta > /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Databases/All_"$s".dna.fasta &
done
wait

for s in A D B FpBSx77 C E F Kestrel
do
makeblastdb -dbtype nucl -in /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Databases/All_"$s".dna.fasta -out /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Databases/All_"$s".dna.fasta
done

for s in A D B FpBSx77 C E F Kestrel
do
COMPARATIVE=`echo A D B FpBSx77 C E F Kestrel | perl -pe "s#$s\s+##" | perl -pe "s#\s+$s\b##"`
for c in $COMPARATIVE
do
blastn -outfmt 7 -strand both -evalue 0.001 -perc_identity 97 -num_threads 2 -db /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Databases/All_"$s".dna.fasta -query /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Databases/All_"$c".dna.fasta -out /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/BLAST_Outs/All_"$s"-v-"$c".dna.outfmt7  &
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 14 ]
do
sleep 10s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
done
wait

function recipricate_all {
f=`echo $i | perl -pe 's#_#-#'`
N=`echo $i | perl -pe 's#.dna.fasta##' | perl -pe 's#.*MSTRG##'`
d=`expr $N / 1000`
mkdir -p /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Merged/"$d"
cat /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/"$s"/"$d"/"$f".dna.fasta > /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Merged/"$d"/Merged"$f".dna.fasta
for c in A D B FpBSx77 C F Kestrel
do
AlignGENE=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/BLAST_Outs/All_"$s"-v-"$c".dna.outfmt7 | awk -v qgene=$i '$2==qgene{print $1; exit}'`
MatchGENE=`grep -v "#" /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/BLAST_Outs/All_"$c"-v-"$s".dna.outfmt7 | awk -v qgene=$i -v tgene=$AlignGENE '$2==tgene && $1==qgene{print $2}$1==tgene{exit}'`

if [ $AlignGENE = $MatchGENE ]
then
MatchFILE=`echo $MatchGENE | perl -pe 's#_#-#'`
mNUM=`echo $MatchGENE | perl -pe 's#.dna.fasta##' | perl -pe 's#.*MSTRG##'`
mDir=`expr $mNUM / 1000`
cat /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/"$c"/"$mDir"/"$MatchFILE".dna.fasta >> /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Merged/"$d"/Merged"$f".dna.fasta
fi
done
}

s=E
ListSubjectGENE=`ls /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/"$s"/*/"$s"-*.dna.fasta | perl -pe 's#.*/##g' | perl -pe 's#\.dna\.fasta##g' | perl -pe 's#-#_#'` 
for i in $ListSubjectGENE
do
recipricate_all &
while [ $JOBS -ge 28 ]
do
sleep 1s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
wait

GENES=`ls /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Merged/*/Merged*.dna.fasta | perl -pe 's#.*/Merged##'`
#Align All Genes
export LD_LIBRARY_PATH=/home/jw5478/.conda/envs/JCOX/lib
for i in $GENES
do
N=`echo $i | perl -pe 's#.dna.fasta##' | perl -pe 's#.*MSTRG##'`
d=`expr $N / 1000`
mkdir -p /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/"$d"
/home/jw5478/Executables/Programs/prank-msa/bin/prank -d=/scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Merged/"$d"/Merged"$i" -iterate=2 -F -codon -quiet -o=/scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/"$d"/Aligned_"$i" &
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 28 ]
do
sleep 1s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Merged

function phylip_fry { 
/home/jw5478/Executables/fasta2relaxedPhylip.pl -f /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/"$d"/Aligned_"$i".best.fas -o /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/Modified/"$i"_SampleAlignedORF.dna.paml
cat /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/Modified/"$i"_SampleAlignedORF.dna.paml >> /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/Modified/AllSampleAlignedORF.dna.paml
rm /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/Modified/"$i"_SampleAlignedORF.dna.paml
}

rm /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/Modified/AllSampleAlignedORF.dna.paml
for i in $(ls /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/*/*.dna.fasta.best.fas | perl -pe 's#.*/Aligned_##' | perl -pe 's#\.best\.fas##')
do
N=`echo $i | perl -pe 's#.dna.fasta##' | perl -pe 's#.*MSTRG##'`
d=`expr $N / 1000`
phylip_fry &
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 28 ]
do
sleep 1s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
wait

for d in $(ls /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned | grep -P '[0-9]+$')
do
cd /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned
tar -czf  aligned."$d".tar.gz "$d" &
done
wait
cd /scratch/jw5478/FalconProject/Projects/RefGenome8/Scripts/EvolutionaryPressures/Submissions

for d in $(ls /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned | grep -P '[0-9]+$')
do
rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/"$d" &
done
wait

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs
perl -pe 's#_*R*_*([a-zA-Z0-9]+_MS[A-Z0-9]*)\b#$1#' /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/Modified/AllSampleAlignedORF.dna.paml | perl -pe 's#^A\s+#A            #' | perl -pe 's#^B\s+#B            #' | perl -pe 's#^C\s+#C            #' | perl -pe 's#^D\s+#D            #' | perl -pe 's#^E\s+#E            #' | perl -pe 's#^F\s+#F            #' | perl -pe 's#^FpBSx77\s+#FpBSx77      #' | perl -pe 's#^Kestrel\s+#Kestrel      #' | perl -pe 's#([^^])\s([0-9]+\s[0-9]+)#$1\n$2#' | perl -pe 's#\t#        #g' | awk '/^\s8/{include="yes"}/^\s[0-7]/{include="no"}include=="yes"{print $0}' | grep -P "^[a-zA-Z0-9_]+\b" | perl -pe 's#^[a-zA-Z0-9]+_(MSTRG[0-9]+).*#$1#' | uniq | awk '{print NR"\t"$0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/AllSampleAlignedORF.dna.relabel.paml.key
perl -pe 's#_*R*_*([a-zA-Z0-9]+_MS[A-Z0-9]*)\b#$1#' /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/Modified/AllSampleAlignedORF.dna.paml | perl -pe 's#^A\s+#A            #' | perl -pe 's#^B\s+#B            #' | perl -pe 's#^C\s+#C            #' | perl -pe 's#^D\s+#D            #' | perl -pe 's#^E\s+#E            #' | perl -pe 's#^F\s+#F            #' | perl -pe 's#^FpBSx77\s+#FpBSx77      #' | perl -pe 's#^Kestrel\s+#Kestrel      #' | perl -pe 's#([^^])\s([0-9]+\s[0-9]+)#$1\n$2#' | perl -pe 's#\t#        #g' | awk 'BEGIN{count=0}/^\s8/{include="yes"; count=count+1; print $0; next}/^\s[0-7]/{include="no"}include=="yes"{print $0" "count}' | perl -pe 's#(^[a-zA-Z0-9]+)_MSTRG[0-9]+#$1#' > /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/AllSampleAlignedORF.dna.relabel.paml
cat /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/AllSampleAlignedORF.dna.relabel.paml > /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/AllSampleAlignedORF.dna.relabel.annotate.paml
for g in $(perl -pe 's#\t#__#' /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/AllSampleAlignedORF.dna.relabel.paml.key)
do
GENE_NUM=`echo $g | perl -pe 's#__#\t#' | awk '{print $1}'`
MSTRG=`echo $g | perl -pe 's#__#\t#' | awk '{print $2}'`
sed -i "s# $GENE_NUM\$# $MSTRG#" /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/AllSampleAlignedORF.dna.relabel.annotate.paml
done

cp /scratch/jw5478/FalconProject/Projects/RefGenome8/NucVar/Phylogeny/RAxML/Ulrametric/Ref8Gamma_RAxML_bipartitionsBranchLabels.tre /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/Ref8Gamma_RAxML_bipartitionsBranchLabels.relabel.tre
sed -i 's#Lanner#A#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/Ref8Gamma_RAxML_bipartitionsBranchLabels.relabel.tre
sed -i 's#Barbay#B#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/Ref8Gamma_RAxML_bipartitionsBranchLabels.relabel.tre
sed -i 's#Peregrine#C#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/Ref8Gamma_RAxML_bipartitionsBranchLabels.relabel.tre
sed -i 's#Saker#D#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/Ref8Gamma_RAxML_bipartitionsBranchLabels.relabel.tre
sed -i 's#Gyr-1#E#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/Ref8Gamma_RAxML_bipartitionsBranchLabels.relabel.tre
sed -i 's#Gyr-2#F#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/Ref8Gamma_RAxML_bipartitionsBranchLabels.relabel.tre
sed -i 's#Bl. Shaheen#FpBSx77#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/Ref8Gamma_RAxML_bipartitionsBranchLabels.relabel.tre
sed -i 's#Kestrel#Kestrel#g' /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/Ref8Gamma_RAxML_bipartitionsBranchLabels.relabel.tre

rm -r /scratch/jw5478/Shared/FalconProject/PAML_Inputs
cp -r /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs /scratch/jw5478/Shared/FalconProject/PAML_Inputs
chmod 755 /scratch/jw5478/Shared/FalconProject/PAML_Inputs
chmod 755 /scratch/jw5478/Shared/FalconProject/PAML_Inputs/*

#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --mem=118GB
#SBATCH -p serial
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=28
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

module purge
module load jdk

FalconidaeMito="/scratch/jw5478/FalconProject/Projects/RefGenome8/Mitochondria/NUMT"
HMMER="/home/jw5478/Executables/Programs/HMMER/bin"
Mafft="/home/jw5478/Executables/Programs/Mafft"
DATE=`ls $FalconidaeMito/Downloads/Sequences/FalconidaeMitoGenome_*.orient.fasta | perl -pe 's#.*FalconidaeMitoGenome_(.*)\.orient.fasta#\1#'`
FalcoMito="/scratch/jw5478/FalconProject/Projects/RefGenome8/Mitochondria/Mitogenome"
rm -r $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees
rm -r $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/FinalTrees
mkdir $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees
mkdir $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/FinalTrees

for f in $(grep ">" $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/AllFalco_16KBph2.diploid.numts.centID80.clust.fasta | perl -pe 's#>##' | perl -pe 's#;.*##' | sort -V)
do
OUTGROUP=`awk 'BEGIN{count=0}/>/{count=count+1}/>AF338715.1/{exit}END{print count}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Aligned/"$f"AllFalco_16KBph2.diploid.numts.linsi.fasta | perl -pe 's#([0-9]+)#$1\_NUMT#'`
cd $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/
rm -r $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/"$f"
mkdir $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/"$f"
cd $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/"$f"
awk 'BEGIN{count=0}/^>/{count=count+1; print ">"count"_NUMT"}!/^>/{print $0}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Aligned/"$f"AllFalco_16KBph2.diploid.numts.linsi.fasta > $FalconidaeMito/NUMT_Phylogeny/Alignment/"$f"AllFalco_16KBph2.diploid.numts.linsi.num.fasta

java -cp /home/jw5478/Executables/Programs/readseq.jar run -inform=8 -f=12 -p $FalconidaeMito/NUMT_Phylogeny/Alignment/"$f"AllFalco_16KBph2.diploid.numts.linsi.num.fasta > $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Aligned/"$f"AllFalco_16KBph2.diploid.numts.linsi.num.phylip
cp $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Aligned/"$f"AllFalco_16KBph2.diploid.numts.linsi.num.phylip $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/"$f"/"$f"AllFalco_16KBph2.diploid.numts.linsi.num.phylip

/home/jw5478/Executables/Programs/standard-RAxML/raxmlHPC-PTHREADS -s $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/"$f"/"$f"AllFalco_16KBph2.diploid.numts.linsi.num.phylip -p 31415927 -o $OUTGROUP -n "$f"NUMTsRef8T1 -m GTRGAMMA -# 5 -T 28
/home/jw5478/Executables/Programs/standard-RAxML/raxmlHPC-PTHREADS -s $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/"$f"/"$f"AllFalco_16KBph2.diploid.numts.linsi.num.phylip -p 31415927 -o $OUTGROUP -n "$f"NUMTsRef8T2 -m GTRGAMMA -b 271828 -# 20 -T 28
/home/jw5478/Executables/Programs/standard-RAxML/raxmlHPC-PTHREADS -s $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/"$f"/"$f"AllFalco_16KBph2.diploid.numts.linsi.num.phylip -f b -t RAxML_bestTree."$f"NUMTsRef8T1 -z RAxML_bootstrap."$f"NUMTsRef8T2 -p 31415927 -o $OUTGROUP -n "$f"NUMTsRef8T3 -m GTRGAMMA -T 28
perl -pe 's#(:0\.[0-9]+)\[([0-9]+)+\]#$2$1#g' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/"$f"/RAxML_bipartitionsBranchLabels."$f"NUMTsRef8T3 > $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/"$f"/RAxML_bipartitionsBranchLabels."$f".tre
COUNT=0
SEQUENCE=`grep ">" $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Aligned/"$f"AllFalco_16KBph2.diploid.numts.linsi.fasta | perl -pe 's#>##' | perl -pe 's#\s+-##g' | perl -pe 's#\b\s+\b#_#g' | perl -pe 's#,.*##' | perl -pe 's#_[a-zA-Z0-9]+:.*##'`
for n in $SEQUENCE
do
NAME=`echo $n | perl -pe 's#^A[:_]#Lanner#g' | perl -pe 's#^B[:_]#Barbary#g' | perl -pe 's#^C[:_]#Peregrine#g' | perl -pe 's#^D[:_]#Saker#g' | perl -pe 's#^E[:_]#Gyr-1#g' | perl -pe 's#^F[:_]#Gyr-2#g' | perl -pe 's#^FpBSx77[:_]#Bl.Shaheen#g' | perl -pe 's#^Kestrel[:_]#Kestrel#g'` 
COUNT=`expr $COUNT + 1` 
sed -r -i "s#\b$COUNT\_NUMT#$NAME#g" $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/"$f"/RAxML_bipartitionsBranchLabels."$f".tre
done
cp $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/"$f"/RAxML_bipartitionsBranchLabels."$f".tre $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/FinalTrees/RAxML_bipartitionsBranchLabels."$f".tre
done

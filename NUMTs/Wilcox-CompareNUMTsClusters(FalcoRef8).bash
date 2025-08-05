#!/bin/bash
#SBATCH --time=2:00:00
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

FalconidaeMito="/scratch/jw5478/FalconProject/Projects/RefGenome8/Mitochondria/NUMT"
rm -r $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTrees
rm -r $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables
mkdir $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTrees
mkdir $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables

function cluster_compile {
perl -pe 's#(:[0-9]\.[0-9]+)\[([0-9]+)+\]#$2$1#g' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/"$f"/RAxML_bipartitionsBranchLabels."$f"NUMTsRef8T3 > $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTrees/RAxML_bipartitionsBranchLabels."$f".tre
COUNT=0
NUMBER=`echo $f | perl -pe 's#Cluster##'`
SEQUENCE=`grep ">" $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Aligned/"$f"AllFalco_16KBph2.diploid.numts.linsi.fasta | perl -pe 's#>##' | perl -pe 's#\s+-##g' | perl -pe 's#\b\s+\b#_#g' | perl -pe 's#,.*##' | perl -pe 's#_[a-zA-Z0-9]+:.*##'`
echo -e "NAME\tn\tCOUNT\tgenSCAFFOLD\tgenPOS1\tgenPOS2\tsynSCAFFOLD\tsynPOS1\tsynPOS2\tTYPE\tBASEPAIRS\tSYNT" > $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.txt
for n in $SEQUENCE
do
TYPE=NUMT
NAME=`grep "$n" $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Aligned/"$f"AllFalco_16KBph2.diploid.numts.linsi.fasta | perl -pe 's#>##'| perl -pe 's#^A[:_]#Lanner#g' | perl -pe 's#^B[:_]#Barbary#g' | perl -pe 's#^C[:_]#Peregrine#g' | perl -pe 's#^D[:_]#Saker#g' | perl -pe 's#^E[:_]#Gyr-1#g' | perl -pe 's#^F[:_]#Gyr-2#g' | perl -pe 's#^FpBSx77[:_]#Bl.Shaheen#g' | perl -pe 's#^Kestrel[:_]#Kestrel#g' | perl -pe 's#:#@#'`
genSCAFFOLD=`echo $NAME | perl -pe 's#@.*##' | perl -pe 's#__#\t#g' | awk '{print $1}'`
genPOS1=`echo $NAME | perl -pe 's#@.*##' | perl -pe 's#__#\t#g' | awk '{print $2}'`
genPOS2=`echo $NAME | perl -pe 's#@.*##' | perl -pe 's#__#\t#g' | awk '{print $3}'`
BASEPAIRS=`echo -e "$genPOS1\t$genPOS2" | awk '$2>$1{basepairs=$2-$1}$1>$2{basepairs=$1-$2}END{print basepairs}'`

synSCAFFOLD=NA
synPOS1=NA
synPOS2=NA
ALIGNED=NA
SYN=`echo $NAME | perl -pe 's#.*@##'`
if [ "$SYN" ]
then
synSCAFFOLD=`echo $NAME | perl -pe 's#.*@##' | perl -pe 's#__#\t#g' | awk '{print $1}'`
synPOS1=`echo $NAME | perl -pe 's#.*@##' | perl -pe 's#__#\t#g' | awk '{print $2}'`
synPOS2=`echo $NAME | perl -pe 's#.*@##' | perl -pe 's#__#\t#g' | awk '{print $3}'`
fi

if [ -z "$NAME" ]
then
NAME=`echo "$n" | perl -pe 's#>##'| perl -pe 's#^A[:_]#Lanner#g' | perl -pe 's#^B[:_]#Barbary#g' | perl -pe 's#^C[:_]#Peregrine#g' | perl -pe 's#^D[:_]#Saker#g' | perl -pe 's#^E[:_]#Gyr-1#g' | perl -pe 's#^F[:_]#Gyr-2#g' | perl -pe 's#^FpBSx77[:_]#Bl.Shaheen#g' | perl -pe 's#^Kestrel[:_]#Kestrel#g' | perl -pe 's#:#@#'` 
genSCAFFOLD=NA
genPOS1=NA
genPOS2=NA
BASEPAIRS=NA
TYPE=MITOGENOME
fi

KESTREL=`echo $n | grep "Kestrel" | perl -pe 's#.*#Kestrel#'`
if [[ $KESTREL == Kestrel ]]
then
NAME=`echo "$n" | perl -pe 's#>##'| perl -pe 's#^A[:_]#Lanner#g' | perl -pe 's#^B[:_]#Barbary#g' | perl -pe 's#^C[:_]#Peregrine#g' | perl -pe 's#^D[:_]#Saker#g' | perl -pe 's#^E[:_]#Gyr-1#g' | perl -pe 's#^F[:_]#Gyr-2#g' | perl -pe 's#^FpBSx77[:_]#Bl.Shaheen#g' | perl -pe 's#^Kestrel[:_]#Kestrel#g' | perl -pe 's#:#@#'` 
synSCAFFOLD=$genSCAFFOLD
synPOS1=$genPOS1
synPOS2=$genPOS2
fi
COUNT=`expr $COUNT + 1` 
sed -r -i "s#\b$COUNT\_NUMT#$NAME#g" $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTrees/RAxML_bipartitionsBranchLabels."$f".tre
SYNT=NA
echo -e "$NAME\t$n\t$COUNT\t$genSCAFFOLD\t$genPOS1\t$genPOS2\t$synSCAFFOLD\t$synPOS1\t$synPOS2\t$TYPE\t$BASEPAIRS\t$SYNT" >> $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.txt
done

cat $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.txt > $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.temp1.txt
SCAFFOLDS=`awk '$10=="NUMT"{print $7}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.txt | sort -V | uniq | perl -pe 's#NA\n##'`
for c in $SCAFFOLDS
do
NUMTs=`awk -v scaffold="$c" '$7==scaffold{print $1}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.txt`
synCOUNT=0
for n in $NUMTs
do
synSCAFFOLD="$c"
nSTART=`awk -v numt="$n" '$1==numt{print $8}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.txt`
nEND=`awk -v numt="$n" '$1==numt{print $9}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.txt`
synCOUNT=`awk -v scaffold=$synSCAFFOLD -v pos1=$nSTART -v pos2=$nEND '$7==scaffold && $8>=pos1 && $9<=pos2 && $12!="NA"{synt=$12; print synt}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.txt | perl -pe 's#.*p##' | sort -n | uniq`
if [ -z $synCOUNT ]
then
synCOUNT=`expr $synCOUNT + 1`
fi
SYNT=`echo $synSCAFFOLD\p$synCOUNT`
awk -v scaffold=$synSCAFFOLD -v pos1=$nSTART -v pos2=$nEND -v synt=$SYNT 'BEGIN{OFS=FS="\t"}$7==scaffold && $8>=pos1 && $9<=pos2{$12=synt; print $0; next}{print $0}' $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.temp1.txt > $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.temp2.txt
cat $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.temp2.txt > $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.temp1.txt
done
done
cat $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.temp1.txt > $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.txt
rm $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/Trees/UltimateTables/RAxML_bipartitionsBranchLabels."$f"_table.temp*.txt
}

for f in $(grep ">" $FalconidaeMito/NUMT_Phylogeny/Alignment/Clusters/AllFalco_16KBph2.diploid.numts.centID80.clust.fasta | perl -pe 's#>##' | perl -pe 's#;.*##' | sort -V)
do
cluster_compile &
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 28 ]
do
sleep 3s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
wait

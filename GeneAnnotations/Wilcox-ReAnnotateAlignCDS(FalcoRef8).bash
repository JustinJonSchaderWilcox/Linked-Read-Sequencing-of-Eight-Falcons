#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --mem=118GB
#SBATCH -p serial
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

rm -r /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs
mkdir /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs
perl -pe 's#_*R*_*([a-zA-Z0-9]+_MS[A-Z0-9]*)\b#$1#' /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/Modified/AllSampleAlignedORF.dna.paml | perl -pe 's#^A\s+#A            #' | perl -pe 's#^B\s+#B            #' | perl -pe 's#^C\s+#C            #' | perl -pe 's#^D\s+#D            #' | perl -pe 's#^E\s+#E            #' | perl -pe 's#^F\s+#F            #' | perl -pe 's#^FpBSx77\s+#FpBSx77      #' | perl -pe 's#^Kestrel\s+#Kestrel      #' | perl -pe 's#([^^])\s([0-9]+\s[0-9]+)#$1\n$2#' | perl -pe 's#\t#        #g' | awk '/^\s8/{include="yes"}/^\s[0-7]/{include="no"}include=="yes"{print $0}' | grep -P "^[a-zA-Z0-9_]+\b" | perl -pe 's#^[a-zA-Z0-9]+_(MSTRG[0-9]+).*#$1#' | uniq | awk '{print NR"\t"$0}' > /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/AllSampleAlignedORF.dna.relabel.paml.key
perl -pe 's#_*R*_*([a-zA-Z0-9]+_MS[A-Z0-9]*)\b#$1#' /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/CDS/Aligned/Modified/AllSampleAlignedORF.dna.paml | perl -pe 's#^A\s+#A            #' | perl -pe 's#^B\s+#B            #' | perl -pe 's#^C\s+#C            #' | perl -pe 's#^D\s+#D            #' | perl -pe 's#^E\s+#E            #' | perl -pe 's#^F\s+#F            #' | perl -pe 's#^FpBSx77\s+#FpBSx77      #' | perl -pe 's#^Kestrel\s+#Kestrel      #' | perl -pe 's#([^^])\s([0-9]+\s[0-9]+)#$1\n$2#' | perl -pe 's#\t#        #g' | awk 'BEGIN{count=0}/^\s8/{include="yes"; count=count+1; print $0; next}/^\s[0-7]/{include="no"}include=="yes"{print $0" "count}' | perl -pe 's#(^[a-zA-Z0-9]+)_MSTRG[0-9]+#$1#' > /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/AllSampleAlignedORF.dna.relabel.paml
cat /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/AllSampleAlignedORF.dna.relabel.paml > /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/AllSampleAlignedORF.dna.relabel.annotate.paml
for g in $(perl -pe 's#\t#__#' /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/AllSampleAlignedORF.dna.relabel.paml.key)
do
GENE_NUM=`echo $g | perl -pe 's#__#\t#' | awk '{print $1}'`
MSTRG=`echo $g | perl -pe 's#__#\t#' | awk '{print $2}'`
sed -E -i "s#([^0-9]) $GENE_NUM\$#\1 $MSTRG#" /scratch/jw5478/FalconProject/Projects/RefGenome8/EvolutionaryPressures/PAML_Inputs/AllSampleAlignedORF.dna.relabel.annotate.paml
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

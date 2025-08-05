sink(file = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Notebooks/UltimateFigurePrepOutput.txt", append = FALSE, type = c("output", "message"),split = FALSE)
#Make to Compare Variants Methods
#Import Data
setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis/SNP_Analysis/Final_Analyses/SNP_Counts")
SelfAlignVariantCounts=read.delim("MethodsCompVariantSummary.txt", row.names=1, header=TRUE, sep="\t")

#Calculate Shared and unique variants
SelfAlignVariantCounts$MummerSNVs_Shared=SelfAlignVariantCounts$mumHET_SNVs - SelfAlignVariantCounts$UmumHET_SNVs
SelfAlignVariantCounts$MummerIndels_Shared=SelfAlignVariantCounts$mumHET_INDELs - SelfAlignVariantCounts$UmumHET_INDELs

SelfAlignVariantCounts$lrSNV_Shared=SelfAlignVariantCounts$lrHET_SNVs - SelfAlignVariantCounts$UlrHET_SNVs
SelfAlignVariantCounts$lrInDel_Shared=SelfAlignVariantCounts$lrHET_INDELs - SelfAlignVariantCounts$UlrHET_INDELs

#Change Labels
SelfAlignVariantCounts$ID=row.names(SelfAlignVariantCounts)
SelfAlignVariantCounts$ID[SelfAlignVariantCounts$ID=="A"]="Lanner"
SelfAlignVariantCounts$ID[SelfAlignVariantCounts$ID=="B"]="Barbary"
SelfAlignVariantCounts$ID[SelfAlignVariantCounts$ID=="C"]="Peregrine"
SelfAlignVariantCounts$ID[SelfAlignVariantCounts$ID=="D"]="Saker"
SelfAlignVariantCounts$ID[SelfAlignVariantCounts$ID=="E"]="Gyr-1"
SelfAlignVariantCounts$ID[SelfAlignVariantCounts$ID=="F"]="Gyr-2"
SelfAlignVariantCounts$ID[SelfAlignVariantCounts$ID=="FpBSx77"]="Bl. Shaheen"
SelfAlignVariantCounts$ID[SelfAlignVariantCounts$ID=="Kestrel"]="C. Kestrel"

drops=c("mumHET_SNVs","mumHET_INDELs","lrHET_SNVs","lrHET_INDELs")
SelfAlignVariantCounts=SelfAlignVariantCounts[ , !(names(SelfAlignVariantCounts) %in% drops)]

library(reshape2)
library(ggplot2)
library(RColorBrewer)

SelfAlignVariantCounts.m=melt(SelfAlignVariantCounts, id.vars="ID")

#Category names for clusters
SelfAlignVariantCounts.m$cat = as.character(SelfAlignVariantCounts.m$variable)

SelfAlignVariantCounts.m$cat[SelfAlignVariantCounts.m$cat =='MummerSNVs_Shared']="Self Mummer"
SelfAlignVariantCounts.m$cat[SelfAlignVariantCounts.m$cat =='UmumHET_SNVs']="Self Mummer"
SelfAlignVariantCounts.m$cat[SelfAlignVariantCounts.m$cat =='MummerIndels_Shared']="Self Mummer"
SelfAlignVariantCounts.m$cat[SelfAlignVariantCounts.m$cat =='UmumHET_INDELs']="Self Mummer"

SelfAlignVariantCounts.m$cat[SelfAlignVariantCounts.m$cat =='lrSNV_Shared']="Self Longranger"
SelfAlignVariantCounts.m$cat[SelfAlignVariantCounts.m$cat =='UlrHET_SNVs']="Self Longranger"
SelfAlignVariantCounts.m$cat[SelfAlignVariantCounts.m$cat =='lrInDel_Shared']="Self Longranger"
SelfAlignVariantCounts.m$cat[SelfAlignVariantCounts.m$cat =='UlrHET_INDELs']="Self Longranger"

SelfAlignVariantCounts.m$cat[SelfAlignVariantCounts.m$cat =='KEST_SNVsHet']="Kestrel Mummer"
SelfAlignVariantCounts.m$cat[SelfAlignVariantCounts.m$cat =='KEST_INDELs']="Kestrel Mummer"

#Type Names for fill
SelfAlignVariantCounts.m$variable=as.character(SelfAlignVariantCounts.m$variable)
SelfAlignVariantCounts.m$variable[SelfAlignVariantCounts.m$variable =='MummerSNVs_Shared']="SNVs, Shared"
SelfAlignVariantCounts.m$variable[SelfAlignVariantCounts.m$variable =='UmumHET_SNVs']="SNVs, Unique"
SelfAlignVariantCounts.m$variable[SelfAlignVariantCounts.m$variable =='MummerIndels_Shared']="Indels, Shared"
SelfAlignVariantCounts.m$variable[SelfAlignVariantCounts.m$variable =='UmumHET_INDELs']="Indels, Unique"

SelfAlignVariantCounts.m$variable[SelfAlignVariantCounts.m$variable =='lrSNV_Shared']="SNVs, Shared"
SelfAlignVariantCounts.m$variable[SelfAlignVariantCounts.m$variable =='UlrHET_SNVs']="SNVs, Unique"
SelfAlignVariantCounts.m$variable[SelfAlignVariantCounts.m$variable =='lrInDel_Shared']="Indels, Shared"
SelfAlignVariantCounts.m$variable[SelfAlignVariantCounts.m$variable =='UlrHET_INDELs']="Indels, Unique"

SelfAlignVariantCounts.m$variable[SelfAlignVariantCounts.m$variable =='KEST_SNVsHet']="SNVs, Kestrel"
SelfAlignVariantCounts.m$variable[SelfAlignVariantCounts.m$variable =='KEST_INDELs']="Indels, Kestrel"

#Make cluster and stacked bar charts

mycolors = c("#202020", "#808080", "#B0B0B0", "#00008b", "#9400d3", "#ff8c00")


SelfAlignVariantCounts.m$variable=factor(SelfAlignVariantCounts.m$variable, levels = c("Indels, Shared", "Indels, Unique", "Indels, Kestrel", "SNVs, Kestrel", "SNVs, Shared", "SNVs, Unique"))
SelfAlignVariantCounts.m$ID=factor(SelfAlignVariantCounts.m$ID, levels = c("Barbary", "Bl. Shaheen", "Peregrine", "Lanner", "Gyr-1", "Gyr-2", "Saker", "C. Kestrel"))

ggplot(SelfAlignVariantCounts.m, aes(cat, value, fill = variable)) + geom_bar(stat="identity") + scale_fill_manual(values=mycolors) + theme(text= element_text(size=12)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, color = "black"), legend.text=element_text(size=6), legend.title = element_blank()) + labs(x="",y="Variants") + guides(color = guide_legend(override.aes = list(size = 0.5))) + facet_wrap(vars(ID), nrow=2, ncol=4) + theme(strip.text = element_text(face = "bold"))
ggsave("suppFig2_VariantMethodsCompare.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Supplementary", scale = 1, width = 7, height = 6.5, units = c("in"), dpi =1200, limitsize = TRUE)
###########################################################
#Make Figure 1 and Unique Structural Variant Counts
#Read in File
setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis/RepeatAnalysis/VariantSummaries")
UniqueVariantCounts=read.delim("strVariantSummary.txt", header=TRUE, row.names=1, sep="\t")
#Stacked Bar Chart
##Scale to columns to 100% and reorder columns
#Rename
UniqueVariantCounts$Genome=as.character(row.names(UniqueVariantCounts))
UniqueVariantCounts$Genome[UniqueVariantCounts$Genome=="A"]="Lanner"
UniqueVariantCounts$Genome[UniqueVariantCounts$Genome=="B"]="Barbary"
UniqueVariantCounts$Genome[UniqueVariantCounts$Genome=="C"]="Peregrine"
UniqueVariantCounts$Genome[UniqueVariantCounts$Genome=="D"]="Saker"
UniqueVariantCounts$Genome[UniqueVariantCounts$Genome=="E"]="Gyr-1"
UniqueVariantCounts$Genome[UniqueVariantCounts$Genome=="F"]="Gyr-2"
UniqueVariantCounts$Genome[UniqueVariantCounts$Genome=="FpBSx77"]="Bl. Shaheen"
UniqueVariantCounts$Genome[UniqueVariantCounts$Genome=="Kestrel"]="C. Kestrel"

library(reshape2)
library(ggplot2)
library(RColorBrewer)
UniqueVariantCounts.m=melt(UniqueVariantCounts, id.vars="Genome")
UniqueVariantCounts.m=UniqueVariantCounts.m[grepl("COUNT", UniqueVariantCounts.m$variable), ]

UniqueVariantCounts.m$variable=sub('^DEL_COUNT', "Deletions", UniqueVariantCounts.m$variable)
UniqueVariantCounts.m$variable=sub('^INS_COUNT', "Insertions", UniqueVariantCounts.m$variable)
UniqueVariantCounts.m$variable=sub('^INV_COUNT', "Inversions", UniqueVariantCounts.m$variable)
UniqueVariantCounts.m$variable=sub('^SegDEL_COUNT', "Segmental Deletions", UniqueVariantCounts.m$variable)
UniqueVariantCounts.m$variable=sub('^SegDUP_COUNT', "Segmental Duplications", UniqueVariantCounts.m$variable)

UniqueVariantCounts.m$Genome=factor(UniqueVariantCounts.m$Genome, levels = c("Barbary", "Bl. Shaheen", "Peregrine", "Lanner", "Gyr-1", "Gyr-2", "Saker", "C. Kestrel"))
mycolors = c("#404040", "#C0C0C0", "#48BF91", "#508030", "#70C040", "#90F050", "#FEE391", "#FE9929", "#B15928", "#662506", "#D4B9DA", "#DF65B0", "#CE1256", "#67001F")
ggplot(UniqueVariantCounts.m, aes(Genome, value, fill = variable)) + geom_bar(stat="identity") + scale_fill_manual(values=mycolors) + theme(text= element_text(size=12)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, color = "black", face="bold"), legend.text=element_text(size=10), legend.title = element_blank()) + labs(x="",y="Unique Variant Counts") + guides(color = guide_legend(override.aes = list(size = 2)))
ggsave("legendFig1A_UniqueSV.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Intermediate", scale = 1, width = 7, height = 5.0, units = c("in"), dpi =1200, limitsize = TRUE)

ggplot(UniqueVariantCounts.m, aes(Genome, value, fill = variable)) + geom_bar(stat="identity") + scale_fill_manual(values=mycolors) + theme(text= element_text(size=12)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, color = "black", face="bold"), legend.text=element_text(size=6), legend.title = element_blank(), legend.position="none") + labs(x="",y="Unique Variant Counts") + guides(color=FALSE)
ggsave("Fig1A_UniqueSV.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Intermediate", scale = 1, width = 7, height = 5.0, units = c("in"), dpi =1200, limitsize = TRUE)

SV_Count.anova=aov(log2(value+1)~variable + Genome, UniqueVariantCounts.m)
SV_Count_Tukey=TukeyHSD(x=SV_Count.anova, "variable", conf.level=0.95)
summary(SV_Count.anova)
SV_Count_Tukey

#Read in File
setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis/RepeatAnalysis/VariantSummaries")
UniqueVariantBP=read.delim("strVariantSummary.txt", header=TRUE, row.names=1, sep="\t")
#Stacked Bar Chart
##Scale to columns to 100% and reorder columns
#Rename
UniqueVariantBP$Genome=as.character(row.names(UniqueVariantBP))
UniqueVariantBP$Genome[UniqueVariantBP$Genome=="A"]="Lanner"
UniqueVariantBP$Genome[UniqueVariantBP$Genome=="B"]="Barbary"
UniqueVariantBP$Genome[UniqueVariantBP$Genome=="C"]="Peregrine"
UniqueVariantBP$Genome[UniqueVariantBP$Genome=="D"]="Saker"
UniqueVariantBP$Genome[UniqueVariantBP$Genome=="E"]="Gyr-1"
UniqueVariantBP$Genome[UniqueVariantBP$Genome=="F"]="Gyr-2"
UniqueVariantBP$Genome[UniqueVariantBP$Genome=="FpBSx77"]="Bl. Shaheen"
UniqueVariantBP$Genome[UniqueVariantBP$Genome=="Kestrel"]="C. Kestrel"

library(reshape2)
library(ggplot2)
library(RColorBrewer)
UniqueVariantBP.m=melt(UniqueVariantBP, id.vars="Genome")
UniqueVariantBP.m=UniqueVariantBP.m[grepl("_BP", UniqueVariantBP.m$variable), ]

UniqueVariantBP.m$variable=sub('^DEL_BP', "Deletions", UniqueVariantBP.m$variable)
UniqueVariantBP.m$variable=sub('^INS_BP', "Insertions", UniqueVariantBP.m$variable)
UniqueVariantBP.m$variable=sub('^INV_BP', "Inversions", UniqueVariantBP.m$variable)
UniqueVariantBP.m$variable=sub('^SegDEL_BP', "Segmental Deletions", UniqueVariantBP.m$variable)
UniqueVariantBP.m$variable=sub('^SegDUP_BP', "Segmental Duplications", UniqueVariantBP.m$variable)

UniqueVariantBP.m$variable=factor(UniqueVariantBP.m$variable, levels=unique(UniqueVariantBP.m$variable))
UniqueVariantBP.m$Genome=factor(UniqueVariantBP.m$Genome, levels = c("Barbary", "Bl. Shaheen", "Peregrine", "Lanner", "Gyr-1", "Gyr-2", "Saker", "C. Kestrel"))
mycolors = c("#404040", "#C0C0C0", "#48BF91", "#508030", "#70C040", "#90F050", "#FEE391", "#FE9929", "#B15928", "#662506", "#D4B9DA", "#DF65B0", "#CE1256", "#67001F")
ggplot(UniqueVariantBP.m, aes(Genome, log10(value), fill = variable)) + geom_bar(stat="identity") + scale_fill_manual(values=mycolors) + theme(text= element_text(size=12)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, color = "black", face="bold"), legend.text=element_text(size=6), legend.title = element_blank(), legend.position="none") + labs(x="",y="log(Unique Variant Base Pairs)") + guides(color=FALSE)
ggsave("Fig1B_UniqueSVbp.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Intermediate", scale = 1, width = 7, height = 5.0, units = c("in"), dpi =1200, limitsize = TRUE)

UniqueVariant_AvgSize=data.frame(Genome=UniqueVariantBP$Genome,DEL_AVG=UniqueVariantBP$DEL_BP/UniqueVariantBP$DEL_COUNT,INS_AVG=UniqueVariantBP$INS_BP/UniqueVariantBP$INS_COUNT,INV_AVG=UniqueVariantBP$INV_BP/UniqueVariantBP$INV_COUNT,SegDEL_AVG=UniqueVariantBP$SegDEL_BP/UniqueVariantBP$SegDEL_COUNT,SegDUP_AVG=UniqueVariantBP$SegDUP_BP/UniqueVariantBP$SegDUP_COUNT,row.names=row.names(UniqueVariantBP))
UniqueVariant_AvgSize.m=melt(UniqueVariant_AvgSize, id.vars="Genome")
UniqueVariant_AvgSize.m$variable=sub('^DEL_AVG', "Deletions", UniqueVariant_AvgSize.m$variable)
UniqueVariant_AvgSize.m$variable=sub('^INS_AVG', "Insertions", UniqueVariant_AvgSize.m$variable)
UniqueVariant_AvgSize.m$variable=sub('^INV_AVG', "Inversions", UniqueVariant_AvgSize.m$variable)
UniqueVariant_AvgSize.m$variable=sub('^SegDEL_AVG', "Segmental Deletions", UniqueVariant_AvgSize.m$variable)
UniqueVariant_AvgSize.m$variable=sub('^SegDUP_AVG', "Segmental Duplications", UniqueVariant_AvgSize.m$variable)

UniqueVariantCounts.m$Type="Count"
UniqueVariantBP.m$Type="Base Pairs"
UniqueVariantBP.m$value=log10(UniqueVariantBP.m$value)
UniqueVariantGraph=rbind(UniqueVariantCounts.m, UniqueVariantBP.m)
UniqueVariantGraph$Type=factor(UniqueVariantGraph$Type, levels = c( "Count", "Base Pairs"))


mycolors = c("#404040", "#C0C0C0", "#48BF91", "#508030", "#70C040", "#90F050", "#FEE391", "#FE9929", "#B15928", "#662506", "#D4B9DA", "#DF65B0", "#CE1256", "#67001F")

ggplot(UniqueVariantGraph, aes(Genome, value, fill = variable)) + geom_bar(stat="identity") + scale_fill_manual(values=mycolors) + theme(text= element_text(size=12)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, color = "black", face="bold"), legend.text=element_text(size=6), legend.title = element_blank(), legend.position="none") + labs(x="",y="") + facet_wrap(vars(Type), nrow=2, scales="free_y") + theme(strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Fig1_UniqueSV.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Intermediate", scale = 1, width = 5, height = 5.0, units = c("in"), dpi =1200, limitsize = TRUE)


#Comparing size differences in unique structural variants
SV_AvgSize.anova=aov(log2(value+1)~variable + Genome, UniqueVariant_AvgSize.m)
SV_AvgSize_Tukey=TukeyHSD(x=SV_AvgSize.anova, "variable", conf.level=0.95)
summary(SV_AvgSize.anova)
SV_AvgSize_Tukey
#without log transformations
SV_AvgSize_notrans.anova=aov(value~variable + Genome, UniqueVariant_AvgSize.m)
SV_AvgSize_notrans_Tukey=TukeyHSD(x=SV_AvgSize_notrans.anova, "variable", conf.level=0.95)
summary(SV_AvgSize_notrans.anova)
SV_AvgSize_notrans_Tukey
#Show min and max genome
UniqueVariantBP$Genome[UniqueVariantBP$SegDEL_BP==min(UniqueVariantBP$SegDEL_BP)]
min(UniqueVariantBP$SegDEL_BP)
UniqueVariantBP$Genome[UniqueVariantBP$SegDEL_BP==max(UniqueVariantBP$SegDEL_BP)]
max(UniqueVariantBP$SegDEL_BP)

UniqueVariantBP$Genome[UniqueVariantBP$INV_BP==min(UniqueVariantBP$INV_BP)]
min(UniqueVariantBP$INV_BP)
UniqueVariantBP$Genome[UniqueVariantBP$INV_BP==max(UniqueVariantBP$INV_BP)]
max(UniqueVariantBP$INV_BP)

UniqueVariantBP$Genome[UniqueVariantBP$INS_BP==min(UniqueVariantBP$INS_BP)]
min(UniqueVariantBP$INS_BP)
UniqueVariantBP$Genome[UniqueVariantBP$INS_BP==max(UniqueVariantBP$INS_BP)]
max(UniqueVariantBP$INS_BP)
#########################################################
##-------Rescale the ith iteration result of PSMC, and make ready for plotting
# file: result file from PSMC
# i.iteration: the ith iteration
# mu: mutation rate
# s: bin size

setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis/SNP_Analysis/Final_Analyses/Final_PSMC/wSexChromosomes/")

psmc.result<-function(file,i.iteration=100,s=100)
{
	X<-scan(file=file,what="",sep="\n",quiet=TRUE)
	
	START<-grep("^RD",X)
	END<-grep("^//",X)
	
	X<-X[START[i.iteration+1]:END[i.iteration+1]]
	
	TR<-grep("^TR",X,value=TRUE)
	RS<-grep("^RS",X,value=TRUE)
	
	write(TR,"temp.psmc.result")
	theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
	
	write(RS,"temp.psmc.result")
	a<-read.table("temp.psmc.result")
	Generation<-as.numeric(a[,3]*theta0/s)
	theta_k<-as.numeric(a[,4]*theta0/s)
	
	file.remove("temp.psmc.result")
	
	n.points<-length(theta_k)
	Divergence<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
		Generation[n.points])
	Pop_Size<-c(as.numeric(rbind(theta_k[-n.points],theta_k[-n.points])),
		theta_k[n.points])
	
	data.frame(Divergence,Pop_Size)
}


refLanner_PSMC=psmc.result("A_16KBph2.diploid.psmc")
refBarbary_PSMC=psmc.result("B_16KBph2.diploid.psmc")
refPeregrine_PSMC=psmc.result("C_16KBph2.diploid.psmc")
refSaker_PSMC=psmc.result("D_16KBph2.diploid.psmc")
refGyrCanadian_PSMC=psmc.result("E_16KBph2.diploid.psmc")
refGyrWrsan_PSMC=psmc.result("F_16KBph2.diploid.psmc")
refBlackShaheen_PSMC=psmc.result("FpBSx77_16KBph2.diploid.psmc")
refKestrel_PSMC=psmc.result("Kestrel_16KBph2.diploid.psmc")

refLanner_PSMC$Boot_Strap=0
refBarbary_PSMC$Boot_Strap=0
refPeregrine_PSMC$Boot_Strap=0
refSaker_PSMC$Boot_Strap=0
refGyrCanadian_PSMC$Boot_Strap=0
refGyrWrsan_PSMC$Boot_Strap=0
refBlackShaheen_PSMC$Boot_Strap=0
refKestrel_PSMC$Boot_Strap=0

refLanner_PSMC=refLanner_PSMC[c(TRUE,FALSE),]
refBarbary_PSMC=refBarbary_PSMC[c(TRUE,FALSE),]
refPeregrine_PSMC=refPeregrine_PSMC[c(TRUE,FALSE),]
refSaker_PSMC=refSaker_PSMC[c(TRUE,FALSE),]
refGyrCanadian_PSMC=refGyrCanadian_PSMC[c(TRUE,FALSE),]
refGyrWrsan_PSMC=refGyrWrsan_PSMC[c(TRUE,FALSE),]
refBlackShaheen_PSMC=refBlackShaheen_PSMC[c(TRUE,FALSE),]
refKestrel_PSMC=refKestrel_PSMC[c(TRUE,FALSE),]

psmc.result<-function(file,i.iteration=100,s=100, bs=100)
{
	X<-scan(file=file,what="",sep="\n",quiet=TRUE)
	
	START<-grep("^RD",X)
	END<-grep("^//",X)
	R=vector(mode="character", length=0)
	for (r in seq(1:bs+1)) {	
	R<-append(R, paste(X[START[i.iteration+1*r]:END[i.iteration+1*r]],r))
	}
	X=R
	rm(R)
	TR<-grep("^TR",X,value=TRUE)
	RS<-grep("^RS",X,value=TRUE)
	
	write(TR,"temp.psmc.result")
	theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
	
	write(RS,"temp.psmc.result")
	a<-read.table("temp.psmc.result")
	Generation<-as.numeric(a[,3]*theta0/s)
	theta_k<-as.numeric(a[,4]*theta0/s)
	Replicate=as.numeric(a[,8])
	file.remove("temp.psmc.result")
	
	n.points<-length(theta_k)
	Divergence<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
		Generation[n.points])
	Pop_Size<-c(as.numeric(rbind(theta_k[-n.points],theta_k[-n.points])),
		theta_k[n.points])
	Boot_Strap<-c(as.numeric(rbind(Replicate[-n.points],Replicate[-n.points])),
		Replicate[n.points])
	data.frame(Divergence,Pop_Size, Boot_Strap)
}

bsLanner_PSMC=psmc.result("A_16KBph2.diploid.combined.psmc")
bsBarbary_PSMC=psmc.result("B_16KBph2.diploid.combined.psmc")
bsPeregrine_PSMC=psmc.result("C_16KBph2.diploid.combined.psmc")
bsSaker_PSMC=psmc.result("D_16KBph2.diploid.combined.psmc")
bsGyrCanadian_PSMC=psmc.result("E_16KBph2.diploid.combined.psmc")
bsGyrWrsan_PSMC=psmc.result("F_16KBph2.diploid.combined.psmc")
bsBlackShaheen_PSMC=psmc.result("FpBSx77_16KBph2.diploid.combined.psmc")
bsKestrel_PSMC=psmc.result("Kestrel_16KBph2.diploid.combined.psmc")

bsLanner_PSMC=bsLanner_PSMC[c(TRUE,FALSE),]
bsBarbary_PSMC=bsBarbary_PSMC[c(TRUE,FALSE),]
bsPeregrine_PSMC=bsPeregrine_PSMC[c(TRUE,FALSE),]
bsSaker_PSMC=bsSaker_PSMC[c(TRUE,FALSE),]
bsGyrCanadian_PSMC=bsGyrCanadian_PSMC[c(TRUE,FALSE),]
bsGyrWrsan_PSMC=bsGyrWrsan_PSMC[c(TRUE,FALSE),]
bsBlackShaheen_PSMC=bsBlackShaheen_PSMC[c(TRUE,FALSE),]
bsKestrel_PSMC=bsKestrel_PSMC[c(TRUE,FALSE),]

Lanner_PSMC=rbind(refLanner_PSMC,bsLanner_PSMC)
Barbary_PSMC=rbind(refBarbary_PSMC,bsBarbary_PSMC)
Peregrine_PSMC=rbind(refPeregrine_PSMC,bsPeregrine_PSMC)
Saker_PSMC=rbind(refSaker_PSMC,bsSaker_PSMC)
GyrCanadian_PSMC=rbind(refGyrCanadian_PSMC,bsGyrCanadian_PSMC)
GyrWrsan_PSMC=rbind(refGyrWrsan_PSMC,bsGyrWrsan_PSMC)
BlackShaheen_PSMC=rbind(refBlackShaheen_PSMC,bsBlackShaheen_PSMC)
Kestrel_PSMC=rbind(refKestrel_PSMC,bsKestrel_PSMC)

#Adding Species
Lanner_PSMC$Species="Lanner"
Barbary_PSMC$Species="Peregrine"
Peregrine_PSMC$Species="Peregrine"
Saker_PSMC$Species="Saker"
GyrCanadian_PSMC$Species="Gyr"
GyrWrsan_PSMC$Species="Gyr"
BlackShaheen_PSMC$Species="Peregrine"
Kestrel_PSMC$Species="Kestrel"

Lanner_PSMC$Genome=Lanner_PSMC$Species
Barbary_PSMC$Genome[Barbary_PSMC$Species=="Peregrine"]="Barbary"
Peregrine_PSMC$Genome=Peregrine_PSMC$Species
Saker_PSMC$Genome=Saker_PSMC$Species
GyrCanadian_PSMC$Genome[GyrCanadian_PSMC$Species=="Gyr"]="Gyr-1"
GyrWrsan_PSMC$Genome[GyrWrsan_PSMC$Species=="Gyr"]="Gyr-2"
BlackShaheen_PSMC$Genome[BlackShaheen_PSMC$Species=="Peregrine"]="Bl. Shaheen"
Kestrel_PSMC$Genome[Kestrel_PSMC$Species=="Kestrel"]="C. Kestrel"
library(ggplot2)

MergedPSMC=as.data.frame(rbind(Lanner_PSMC, Barbary_PSMC, Peregrine_PSMC, Saker_PSMC, GyrCanadian_PSMC, GyrWrsan_PSMC, BlackShaheen_PSMC,Kestrel_PSMC))
MergedPSMC$Genome=factor(MergedPSMC$Genome, levels = c("Barbary", "Bl. Shaheen", "Peregrine", "Lanner", "Gyr-1", "Gyr-2", "Saker", "C. Kestrel"))

MergedPSMC.pruned=MergedPSMC[MergedPSMC$Genome %in% c("Peregrine", "Lanner", "Gyr-1", "Saker", "C. Kestrel"),]
MergedPSMC.pruned=MergedPSMC.pruned[MergedPSMC.pruned$Divergence<=0.001634,]

MergedPSMC.pruned$Genome=factor(MergedPSMC.pruned$Genome, levels = c("Peregrine", "Lanner", "Gyr-1", "Saker", "C. Kestrel"))
MergedPSMC.pruned$gBootstrap=paste(MergedPSMC.pruned$Genome, MergedPSMC.pruned$Boot_Strap)
#Old Pruned MergedPSMC.pruned$gBootstrap=factor(MergedPSMC.pruned$gBootstrap, levels=c(paste(rep("Barbary", 101), seq(100,0)),paste(rep("Bl. Shaheen", 101), seq(100,0)),paste(rep("Peregrine", 101), seq(100,0)),paste(rep("Lanner", 101), seq(100,0)),paste(rep("Gyr-1", 101), seq(100,0)),paste(rep("Gyr-2", 101), seq(100,0)),paste(rep("Kestrel", 101), seq(100,0)),paste(rep("C. Kestrel", 101), seq(100,0))))
MergedPSMC.pruned$gBootstrap=factor(MergedPSMC.pruned$gBootstrap, levels=c(paste(rep("Peregrine", 101), seq(100,0)),paste(rep("Lanner", 101), seq(100,0)),paste(rep("Gyr-1", 101), seq(100,0)),paste(rep("Saker", 101), seq(100,0)),paste(rep("C. Kestrel", 101), seq(100,0))))

#Key to Old Colors: Barbary=Darkblue; BlackShaheen=Green; Peregrine=Firebrick; Lanner=Darkviolet; Gyr-1=Gold4; Gyr-2=DarkOrange; Kestrel=Khaki4; 
#Old mycolors=c(c(rep("#00008b19", 100), "#00008b"), c(rep("#008b0019", 100), "#008b00"),c(rep("#b2222219", 100), "#b22222"), c(rep("#9400d319", 100), "#9400d3"), c(rep("#eec90019", 100), "#eec900"), c(rep("#ff8c0019", 100), "#ff8c00"), c(rep("#8b5a2b19", 100), "#8b5a2b"),c(rep("#8b864e19", 100), "#8b864e"))
mycolors=c(c(rep("#00008b19", 100), "#00008b"), c(rep("#ff8c0019", 100), "#ff8c00"), c(rep("#008b0019", 100), "#008b00"),c(rep("#9400d319", 100), "#9400d3"),c(rep("#b2222219", 100), "#b22222"))
#X-scaled to one million years based on averages across species: species were averaged for one million years divergence, and then an average was taken across species
ggplot(MergedPSMC.pruned) + geom_step(aes(x=Divergence, y=log10(Pop_Size), group=gBootstrap, color=gBootstrap), size=0.8) + scale_color_manual(values=mycolors) + labs(x = "Scaled Pairwise Divergence", y = "log(Scaled Mutation Rate: Theta k)") + theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)), legend.position="none") + guides(color=FALSE)+scale_x_reverse(limits=c(0.001634, -0.0003))
ggsave("Fig4A_PSMC1000k.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Intermediate", scale = 1, width = 7, height = 5.0, units = c("in"), dpi =1200, limitsize = TRUE)
################################################
##-------Rescale the ith iteration result of PSMC, and make ready for plotting
# file: result file from PSMC
# i.iteration: the ith iteration
# mu: mutation rate
# s: bin size
# g: years per generation

setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis/SNP_Analysis/Final_Analyses/Final_PSMC/wSexChromosomes/")

psmc.result<-function(file,i.iteration=100,mu=4.6e-09,s=100,g=3)
{
	X<-scan(file=file,what="",sep="\n",quiet=TRUE)
	
	START<-grep("^RD",X)
	END<-grep("^//",X)
	
	X<-X[START[i.iteration+1]:END[i.iteration+1]]
	
	TR<-grep("^TR",X,value=TRUE)
	RS<-grep("^RS",X,value=TRUE)
	
	write(TR,"temp.psmc.result")
	theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
	N0<-theta0/4/mu/s
	
	write(RS,"temp.psmc.result")
	a<-read.table("temp.psmc.result")
	Generation<-as.numeric(2*N0*a[,3])
	Ne<-as.numeric(N0*a[,4])
	
	file.remove("temp.psmc.result")
	
	n.points<-length(Ne)
	YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
		Generation[n.points])*g
	Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
		Ne[n.points])
	
	data.frame(YearsAgo,Ne)
}

refLanner_PSMC=psmc.result("A_16KBph2.diploid.psmc")
refBarbary_PSMC=psmc.result("B_16KBph2.diploid.psmc")
refPeregrine_PSMC=psmc.result("C_16KBph2.diploid.psmc")
refSaker_PSMC=psmc.result("D_16KBph2.diploid.psmc")
refGyrCanadian_PSMC=psmc.result("E_16KBph2.diploid.psmc")
refGyrWrsan_PSMC=psmc.result("F_16KBph2.diploid.psmc")
refBlackShaheen_PSMC=psmc.result("FpBSx77_16KBph2.diploid.psmc")
refKestrel_PSMC=psmc.result("Kestrel_16KBph2.diploid.psmc")

refLanner_PSMC$Boot_Strap=0
refBarbary_PSMC$Boot_Strap=0
refPeregrine_PSMC$Boot_Strap=0
refSaker_PSMC$Boot_Strap=0
refGyrCanadian_PSMC$Boot_Strap=0
refGyrWrsan_PSMC$Boot_Strap=0
refBlackShaheen_PSMC$Boot_Strap=0
refKestrel_PSMC$Boot_Strap=0

refLanner_PSMC=refLanner_PSMC[c(TRUE,FALSE),]
refBarbary_PSMC=refBarbary_PSMC[c(TRUE,FALSE),]
refPeregrine_PSMC=refPeregrine_PSMC[c(TRUE,FALSE),]
refSaker_PSMC=refSaker_PSMC[c(TRUE,FALSE),]
refGyrCanadian_PSMC=refGyrCanadian_PSMC[c(TRUE,FALSE),]
refGyrWrsan_PSMC=refGyrWrsan_PSMC[c(TRUE,FALSE),]
refBlackShaheen_PSMC=refBlackShaheen_PSMC[c(TRUE,FALSE),]
refKestrel_PSMC=refKestrel_PSMC[c(TRUE,FALSE),]

psmc.result<-function(file,i.iteration=100,mu=4.6e-09,s=100,g=3, bs=100)
{
	X<-scan(file=file,what="",sep="\n",quiet=TRUE)
	
	START<-grep("^RD",X)
	END<-grep("^//",X)
	R=vector(mode="character", length=0)
	for (r in seq(1:bs+1)) {
	
	R<-append(R, paste(X[START[i.iteration+1*r]:END[i.iteration+1*r]],r))
	}
	X=R
	rm(R)
	TR<-grep("^TR",X,value=TRUE)
	RS<-grep("^RS",X,value=TRUE)
	
	write(TR,"temp.psmc.result")
	theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
	N0<-theta0/4/mu/s
	
	write(RS,"temp.psmc.result")
	a<-read.table("temp.psmc.result")
	Generation<-as.numeric(2*N0*a[,3])
	Ne<-as.numeric(N0*a[,4])
	Replicate=as.numeric(a[,8])	
	file.remove("temp.psmc.result")
	
	n.points<-length(Ne)
	YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
		Generation[n.points])*g
	Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
		Ne[n.points])
	Boot_Strap<-c(as.numeric(rbind(Replicate[-n.points],Replicate[-n.points])),
		Replicate[n.points])
	data.frame(YearsAgo,Ne,Boot_Strap)
}

bsLanner_PSMC=psmc.result("A_16KBph2.diploid.combined.psmc")
bsBarbary_PSMC=psmc.result("B_16KBph2.diploid.combined.psmc")
bsPeregrine_PSMC=psmc.result("C_16KBph2.diploid.combined.psmc")
bsSaker_PSMC=psmc.result("D_16KBph2.diploid.combined.psmc")
bsGyrCanadian_PSMC=psmc.result("E_16KBph2.diploid.combined.psmc")
bsGyrWrsan_PSMC=psmc.result("F_16KBph2.diploid.combined.psmc")
bsBlackShaheen_PSMC=psmc.result("FpBSx77_16KBph2.diploid.combined.psmc")
bsKestrel_PSMC=psmc.result("Kestrel_16KBph2.diploid.combined.psmc")

bsLanner_PSMC=bsLanner_PSMC[c(TRUE,FALSE),]
bsBarbary_PSMC=bsBarbary_PSMC[c(TRUE,FALSE),]
bsPeregrine_PSMC=bsPeregrine_PSMC[c(TRUE,FALSE),]
bsSaker_PSMC=bsSaker_PSMC[c(TRUE,FALSE),]
bsGyrCanadian_PSMC=bsGyrCanadian_PSMC[c(TRUE,FALSE),]
bsGyrWrsan_PSMC=bsGyrWrsan_PSMC[c(TRUE,FALSE),]
bsBlackShaheen_PSMC=bsBlackShaheen_PSMC[c(TRUE,FALSE),]
bsKestrel_PSMC=bsKestrel_PSMC[c(TRUE,FALSE),]

Lanner_PSMC=rbind(bsLanner_PSMC,refLanner_PSMC)
Barbary_PSMC=rbind(bsBarbary_PSMC,refBarbary_PSMC)
Peregrine_PSMC=rbind(bsPeregrine_PSMC,refPeregrine_PSMC)
Saker_PSMC=rbind(bsSaker_PSMC,refSaker_PSMC)
GyrCanadian_PSMC=rbind(bsGyrCanadian_PSMC,refGyrCanadian_PSMC)
GyrWrsan_PSMC=rbind(bsGyrWrsan_PSMC,refGyrWrsan_PSMC)
BlackShaheen_PSMC=rbind(bsBlackShaheen_PSMC,refBlackShaheen_PSMC)
Kestrel_PSMC=rbind(bsKestrel_PSMC,refKestrel_PSMC)

#Adding Species
Lanner_PSMC$Species="Lanner"
Barbary_PSMC$Species="Peregrine"
Peregrine_PSMC$Species="Peregrine"
Saker_PSMC$Species="Saker"
GyrCanadian_PSMC$Species="Gyr"
GyrWrsan_PSMC$Species="Gyr"
BlackShaheen_PSMC$Species="Peregrine"
Kestrel_PSMC$Species="Kestrel"

Lanner_PSMC$Genome=Lanner_PSMC$Species
Barbary_PSMC$Genome[Barbary_PSMC$Species=="Peregrine"]="Barbary"
Peregrine_PSMC$Genome=Peregrine_PSMC$Species
Saker_PSMC$Genome=Saker_PSMC$Species
GyrCanadian_PSMC$Genome[GyrCanadian_PSMC$Species=="Gyr"]="Gyr-1"
GyrWrsan_PSMC$Genome[GyrWrsan_PSMC$Species=="Gyr"]="Gyr-2"
BlackShaheen_PSMC$Genome[BlackShaheen_PSMC$Species=="Peregrine"]="Bl. Shaheen"
Kestrel_PSMC$Genome[Kestrel_PSMC$Species=="Kestrel"]="C. Kestrel"
library(ggplot2)

MergedPSMC=as.data.frame(rbind(Lanner_PSMC, Barbary_PSMC, Peregrine_PSMC, Saker_PSMC, GyrCanadian_PSMC, GyrWrsan_PSMC, BlackShaheen_PSMC,Kestrel_PSMC))
MergedPSMC$Genome=factor(MergedPSMC$Genome, levels = c("Barbary", "Bl. Shaheen", "Peregrine", "Lanner", "Gyr-1", "Gyr-2", "Saker", "C. Kestrel"))
MergedPSMC$Boot_Strap=factor(MergedPSMC$Boot_Strap, levels = c(seq(1,100), 0))
MergedPSMC.pruned=MergedPSMC[MergedPSMC$YearsAgo<=600000,]
mycolors=c(rep("grey", 100), "black")
ggplot(MergedPSMC.pruned) + geom_step(aes(x=YearsAgo, y=log10(Ne), group=Boot_Strap, color=Boot_Strap)) + scale_color_manual(values=mycolors) + labs(x = "Years Ago", y = "log(Effective Population Size)") + facet_wrap(vars(Genome), nrow=2, ncol=4, scales="free") + xlim(-100000, 600000) + theme(axis.text.x = element_text(size=rel(0.6)), legend.position="none") + guides(color=FALSE)
ggsave("suppFig5_SuppScaled_PSMC.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Supplementary/", scale = 1, width = 7, height = 5.0, units = c("in"), dpi =1200, limitsize = TRUE)
#######################################
#Make Figure 7, Unique Variant Types
#Read in File
setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis/SNP_Analysis/Final_Analyses/SNP_Counts")
KestrelAlignVariantCounts=read.delim("KestrelAlignVariantCount.txt", header=TRUE, row.names=1, sep="\t")
#Stacked Bar Chart
##Scale to columns to 100% and reorder columns

RelKestrelAlignVariantCounts.scl=as.data.frame(t(scale(t(KestrelAlignVariantCounts), center=FALSE, scale=(colSums(t(KestrelAlignVariantCounts))))))
RelKestrelAlignVariantCounts.scl.ord=RelKestrelAlignVariantCounts.scl[, c("INSERT", "DEL", "AG", "CT", "GA", "TC", "AC", "AT", "CA", "CG", "GC", "GT", "TA", "TG")]
RelKestrelAlignVariantCounts.scl.ord$Genome=row.names(RelKestrelAlignVariantCounts.scl.ord)

#Rename
RelKestrelAlignVariantCounts.scl.ord$Genome[RelKestrelAlignVariantCounts.scl.ord$Genome=="A"]="Lanner"
RelKestrelAlignVariantCounts.scl.ord$Genome[RelKestrelAlignVariantCounts.scl.ord$Genome=="B"]="Barbary"
RelKestrelAlignVariantCounts.scl.ord$Genome[RelKestrelAlignVariantCounts.scl.ord$Genome=="C"]="Peregrine"
RelKestrelAlignVariantCounts.scl.ord$Genome[RelKestrelAlignVariantCounts.scl.ord$Genome=="D"]="Saker"
RelKestrelAlignVariantCounts.scl.ord$Genome[RelKestrelAlignVariantCounts.scl.ord$Genome=="E"]="Gyr-1"
RelKestrelAlignVariantCounts.scl.ord$Genome[RelKestrelAlignVariantCounts.scl.ord$Genome=="F"]="Gyr-2"
RelKestrelAlignVariantCounts.scl.ord$Genome[RelKestrelAlignVariantCounts.scl.ord$Genome=="FpBSx77"]="Bl. Shaheen"
RelKestrelAlignVariantCounts.scl.ord$Genome[RelKestrelAlignVariantCounts.scl.ord$Genome=="Kestrel"]="C. Kestrel"

library(reshape2)
library(ggplot2)
library(RColorBrewer)
RelKestrelAlignVariantCounts.scl.ord.m=melt(RelKestrelAlignVariantCounts.scl.ord, id.vars="Genome")
RelKestrelAlignVariantCounts.scl.ord.m$variable=sub('\\.', " to ", RelKestrelAlignVariantCounts.scl.ord.m$variable)
RelKestrelAlignVariantCounts.scl.ord.m$variable=factor(RelKestrelAlignVariantCounts.scl.ord.m$variable, levels=unique(RelKestrelAlignVariantCounts.scl.ord.m$variable))


library(reshape2)
library(ggplot2)
library(RColorBrewer)
RelKestrelAlignVariantCounts.scl.ord.m=melt(RelKestrelAlignVariantCounts.scl.ord, id.vars="Genome")


RelKestrelAlignVariantCounts.scl.ord.m$variable=sub('^A', "A to ", RelKestrelAlignVariantCounts.scl.ord.m$variable)
RelKestrelAlignVariantCounts.scl.ord.m$variable=sub('^T', "T to ", RelKestrelAlignVariantCounts.scl.ord.m$variable)
RelKestrelAlignVariantCounts.scl.ord.m$variable=sub('^C', "C to ", RelKestrelAlignVariantCounts.scl.ord.m$variable)
RelKestrelAlignVariantCounts.scl.ord.m$variable=sub('^G', "G to ", RelKestrelAlignVariantCounts.scl.ord.m$variable)

RelKestrelAlignVariantCounts.scl.ord.m$variable=sub('INSERT', "Insertion", RelKestrelAlignVariantCounts.scl.ord.m$variable)
RelKestrelAlignVariantCounts.scl.ord.m$variable=sub('DEL', "Deletion", RelKestrelAlignVariantCounts.scl.ord.m$variable)

RelKestrelAlignVariantCounts.scl.ord.m$variable=factor(RelKestrelAlignVariantCounts.scl.ord.m$variable, levels=unique(RelKestrelAlignVariantCounts.scl.ord.m$variable))

#Should group by transition and transvertions and Indel and then color code each group with a different set of colors

#My colors
RelKestrelAlignVariantCounts.scl.ord.m$variable=factor(RelKestrelAlignVariantCounts.scl.ord.m$variable, levels = c("Insertion", "Deletion", "A to G", "T to C", "G to A", "C to T", "C to A", "G to T", "T to G", "A to C", "A to T", "T to A", "C to G", "G to C"))
RelKestrelAlignVariantCounts.scl.ord.m$Genome=factor(RelKestrelAlignVariantCounts.scl.ord.m$Genome, levels = c("Barbary", "Bl. Shaheen", "Peregrine", "Lanner", "Gyr-1", "Gyr-2", "Saker", "C. Kestrel"))

mycolors = c("#404040", "#C0C0C0", "#384828", "#508030", "#70C040", "#90F050", "#FEE391", "#FE9929", "#B15928", "#662506", "#D4B9DA", "#DF65B0", "#CE1256", "#67001F")
ggplot(RelKestrelAlignVariantCounts.scl.ord.m, aes(Genome, value, fill = variable)) + geom_bar(stat="identity") + scale_fill_manual(values=mycolors) + theme(text= element_text(size=12)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, color = "black", face="bold"), legend.text=element_text(size=6), legend.title = element_blank()) + labs(x="",y="Relative Abundance") + guides(color = guide_legend(override.aes = list(size = 0.5)))
ggsave("suppFig6A_RelKestrelAlignVariantCounts.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Intermediate", scale = 1, width = 3.5, height = 5.0, units = c("in"), dpi =1200, limitsize = TRUE)
#######################################
#Assessing Indel Size Differences and Counts
setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis/SNP_Analysis/Final_Analyses/SNP_Counts")

KestrelAlignVariantCounts=read.delim("UniqueKestrelAlignVariantCount.txt", header=TRUE, row.names=1, sep="\t")
HaploidConservedSitesBaseComp=read.delim("HaploidConservedSitesBaseComp.txt", header=TRUE, row.names=1, sep="\t")

KestrelAlignVariantCounts$AT_to_GC=KestrelAlignVariantCounts$AC + KestrelAlignVariantCounts$AG + KestrelAlignVariantCounts$TC + KestrelAlignVariantCounts$TG
KestrelAlignVariantCounts$GC_to_AT=KestrelAlignVariantCounts$GA + KestrelAlignVariantCounts$GT + KestrelAlignVariantCounts$CT + KestrelAlignVariantCounts$CA

KestrelAlignVariantCounts$AT_Total=KestrelAlignVariantCounts$AT + KestrelAlignVariantCounts$AC + KestrelAlignVariantCounts$AG + KestrelAlignVariantCounts$TA + KestrelAlignVariantCounts$TC + KestrelAlignVariantCounts$TG
KestrelAlignVariantCounts$GC_Total=KestrelAlignVariantCounts$GA + KestrelAlignVariantCounts$GT + KestrelAlignVariantCounts$GC + KestrelAlignVariantCounts$CA + KestrelAlignVariantCounts$CT + KestrelAlignVariantCounts$CG

KestrelAlignVariantCounts$AT_Equil=KestrelAlignVariantCounts$AT_Total - KestrelAlignVariantCounts$AT_to_GC + (2*HaploidConservedSitesBaseComp$CombinedAT)
KestrelAlignVariantCounts$GC_Equil=KestrelAlignVariantCounts$GC_Total - KestrelAlignVariantCounts$GC_to_AT + (2*HaploidConservedSitesBaseComp$CombinedGC)

setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis/SNP_Analysis/Final_Analyses/SNP_Counts")
UniqueKestrelAlignIndelSize=read.delim("UniqueKestrelAlignIndelSize.txt", header=TRUE, sep="\t")

UniqueKestrelAlignIndelSize$Genome=as.character(UniqueKestrelAlignIndelSize$Genome)
UniqueKestrelAlignIndelSize$Genome[UniqueKestrelAlignIndelSize$Genome=="A"]="F. biarmicus"
UniqueKestrelAlignIndelSize$Genome[UniqueKestrelAlignIndelSize$Genome=="B"]="F. p. pelegrinoides"
UniqueKestrelAlignIndelSize$Genome[UniqueKestrelAlignIndelSize$Genome=="C"]="F. peregrinus"
UniqueKestrelAlignIndelSize$Genome[UniqueKestrelAlignIndelSize$Genome=="D"]="F. cherrug"
UniqueKestrelAlignIndelSize$Genome[UniqueKestrelAlignIndelSize$Genome=="E"]="F. rusticolus (1)"
UniqueKestrelAlignIndelSize$Genome[UniqueKestrelAlignIndelSize$Genome=="F"]="F. rusticolus (2)"
UniqueKestrelAlignIndelSize$Genome[UniqueKestrelAlignIndelSize$Genome=="FpBSx77"]="F. p. peregrinator"
UniqueKestrelAlignIndelSize$Genome[UniqueKestrelAlignIndelSize$Genome=="Kestrel"]="F. tinnunculus"


#Overall Compared to Kestrel
UniqueKestrelAlignIndelSize$Genome=factor(UniqueKestrelAlignIndelSize$Genome, levels = c("F. tinnunculus", "F. biarmicus", "F. cherrug", "F. p. pelegrinoides", "F. p. peregrinator", "F. peregrinus", "F. rusticolus (1)", "F. rusticolus (2)"))
LM=lm(Length~Type + Genome + Type*Genome + 0, UniqueKestrelAlignIndelSize)
summary(LM)
drop1(LM,.~., test="F")

t.test(UniqueKestrelAlignIndelSize$Length[UniqueKestrelAlignIndelSize$Genome=="F. tinnunculus"]~UniqueKestrelAlignIndelSize$Type[UniqueKestrelAlignIndelSize$Genome=="F. tinnunculus"])


UniqueKestrelAlignIndelSize$Genome=factor(UniqueKestrelAlignIndelSize$Genome, levels = c("F. biarmicus", "F. cherrug", "F. p. pelegrinoides", "F. p. peregrinator", "F. peregrinus", "F. rusticolus (1)", "F. rusticolus (2)", "F. tinnunculus"))


ANOVA=aov(Length~Type + Genome + Type*Genome, UniqueKestrelAlignIndelSize)
TukeyHSD(ANOVA, "Genome")

library(ggplot2)
ggplot(UniqueKestrelAlignIndelSize, aes(x=UniqueKestrelAlignIndelSize$Type, y=UniqueKestrelAlignIndelSize$Length))+geom_boxplot()+xlab("")+labs(y='Length (BP)')+theme(axis.text.x=element_text(angle =- 90, vjust = 0.5, hjust=0, size=12), axis.text.y=element_text(size=12), axis.title=element_text(size=12, face="plain")) + facet_wrap(vars(Genome), nrow=2, ncol=4) + theme(strip.text = element_text(face = "italic"))
ggsave("SuppFig7_IndelSizeDistribution.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Supplementary", scale = 1, width = 7, height = 5.0, units = c("in"), dpi =1200, limitsize = TRUE)

UniqueKestrelAlignVariantCount=read.delim("UniqueKestrelAlignVariantCount.txt", header=TRUE, sep="\t")
UniqueKestrelAlignVariantCount$Genome=as.character(UniqueKestrelAlignVariantCount$X)
UniqueKestrelAlignVariantCount$Genome[UniqueKestrelAlignVariantCount$Genome=="A"]="F. biarmicus"
UniqueKestrelAlignVariantCount$Genome[UniqueKestrelAlignVariantCount$Genome=="B"]="F. p. pelegrinoides"
UniqueKestrelAlignVariantCount$Genome[UniqueKestrelAlignVariantCount$Genome=="C"]="F. peregrinus"
UniqueKestrelAlignVariantCount$Genome[UniqueKestrelAlignVariantCount$Genome=="D"]="F. cherrug"
UniqueKestrelAlignVariantCount$Genome[UniqueKestrelAlignVariantCount$Genome=="E"]="F. rusticolus (1)"
UniqueKestrelAlignVariantCount$Genome[UniqueKestrelAlignVariantCount$Genome=="F"]="F. rusticolus (2)"
UniqueKestrelAlignVariantCount$Genome[UniqueKestrelAlignVariantCount$Genome=="FpBSx77"]="F. p. peregrinator"
UniqueKestrelAlignVariantCount$Genome[UniqueKestrelAlignVariantCount$Genome=="Kestrel"]="F. tinnunculus"

UniqueKestrelAlignVariantCount$Size=UniqueKestrelAlignVariantCount$INSERT + UniqueKestrelAlignVariantCount$DEL

GENOME="F. biarmicus"
print(GENOME)
pbinom(UniqueKestrelAlignVariantCount$DEL[UniqueKestrelAlignVariantCount$Genome==GENOME], UniqueKestrelAlignVariantCount$Size[UniqueKestrelAlignVariantCount$Genome==GENOME], 0.5)

GENOME="F. cherrug"
print(GENOME)
pbinom(UniqueKestrelAlignVariantCount$DEL[UniqueKestrelAlignVariantCount$Genome==GENOME], UniqueKestrelAlignVariantCount$Size[UniqueKestrelAlignVariantCount$Genome==GENOME], 0.5)

GENOME="F. p. pelegrinoides"
print(GENOME)
pbinom(UniqueKestrelAlignVariantCount$DEL[UniqueKestrelAlignVariantCount$Genome==GENOME], UniqueKestrelAlignVariantCount$Size[UniqueKestrelAlignVariantCount$Genome==GENOME], 0.5)

GENOME="F. p. peregrinator"
print(GENOME)
pbinom(UniqueKestrelAlignVariantCount$DEL[UniqueKestrelAlignVariantCount$Genome==GENOME], UniqueKestrelAlignVariantCount$Size[UniqueKestrelAlignVariantCount$Genome==GENOME], 0.5)

GENOME="F. peregrinus"
print(GENOME)
pbinom(UniqueKestrelAlignVariantCount$DEL[UniqueKestrelAlignVariantCount$Genome==GENOME], UniqueKestrelAlignVariantCount$Size[UniqueKestrelAlignVariantCount$Genome==GENOME], 0.5)

GENOME="F. rusticolus (1)"
print(GENOME)
pbinom(UniqueKestrelAlignVariantCount$DEL[UniqueKestrelAlignVariantCount$Genome==GENOME], UniqueKestrelAlignVariantCount$Size[UniqueKestrelAlignVariantCount$Genome==GENOME], 0.5)

GENOME="F. rusticolus (2)"
print(GENOME)
pbinom(UniqueKestrelAlignVariantCount$DEL[UniqueKestrelAlignVariantCount$Genome==GENOME], UniqueKestrelAlignVariantCount$Size[UniqueKestrelAlignVariantCount$Genome==GENOME], 0.5)

GENOME="F. tinnunculus"
print(GENOME)
pbinom(UniqueKestrelAlignVariantCount$INSERT[UniqueKestrelAlignVariantCount$Genome==GENOME], UniqueKestrelAlignVariantCount$Size[UniqueKestrelAlignVariantCount$Genome==GENOME], 0.5)
#######################################
#Rerunning AT-GC Equllibrium with Fixed Sites
#Make Chromosome Annotation Figures
setwd("~/NYU_AD/FalconProject/Manuscripts/FalconRef8/GBE_Revision/UltRevIsochore")
A_revIsochores100KB=read.delim("A_ultChromAnnotatedWindowsSummary.txt", stringsAsFactors = FALSE, header=TRUE, sep="\t")
B_revIsochores100KB=read.delim("B_ultChromAnnotatedWindowsSummary.txt", stringsAsFactors = FALSE, header=TRUE, sep="\t")
C_revIsochores100KB=read.delim("C_ultChromAnnotatedWindowsSummary.txt", stringsAsFactors = FALSE, header=TRUE, sep="\t")
D_revIsochores100KB=read.delim("D_ultChromAnnotatedWindowsSummary.txt", stringsAsFactors = FALSE, header=TRUE, sep="\t")
E_revIsochores100KB=read.delim("E_ultChromAnnotatedWindowsSummary.txt", stringsAsFactors = FALSE, header=TRUE, sep="\t")
F_revIsochores100KB=read.delim("F_ultChromAnnotatedWindowsSummary.txt", stringsAsFactors = FALSE, header=TRUE, sep="\t")
FpBSx77_revIsochores100KB=read.delim("FpBSx77_ultChromAnnotatedWindowsSummary.txt", stringsAsFactors = FALSE, header=TRUE, sep="\t")
Kestrel_revIsochores100KB=read.delim("Kestrel_ultChromAnnotatedWindowsSummary.txt", stringsAsFactors = FALSE, header=TRUE, sep="\t")

F_revIsochores100KB$s="F"

AllrevIsochores100KB=rbind(A_revIsochores100KB, B_revIsochores100KB, C_revIsochores100KB, D_revIsochores100KB, E_revIsochores100KB, F_revIsochores100KB, FpBSx77_revIsochores100KB, Kestrel_revIsochores100KB)
AllrevIsochores100KB.ord=AllrevIsochores100KB[order(AllrevIsochores100KB$s, AllrevIsochores100KB$p),]


rm(A_revIsochores100KB)
rm(B_revIsochores100KB)
rm(C_revIsochores100KB)
rm(D_revIsochores100KB)
rm(E_revIsochores100KB)
rm(F_revIsochores100KB)
rm(FpBSx77_revIsochores100KB)
rm(Kestrel_revIsochores100KB)

Fixed_wCpG=AllrevIsochores100KB.ord[,45:56]
Fixed_nCpG=AllrevIsochores100KB.ord[,77:88]
h_wCpG=AllrevIsochores100KB.ord[,61:72]
h_nCpG=AllrevIsochores100KB.ord[,89:100]

gFixed_wCpG=AllrevIsochores100KB.ord[,c(1,45:56)]
gFixed_nCpG=AllrevIsochores100KB.ord[,c(1,77:88)]
gh_wCpG=AllrevIsochores100KB.ord[,c(1,61:72)]
gh_nCpG=AllrevIsochores100KB.ord[,c(1,89:100)]

#Equillibria for fixed differences with CpG
f_KestrelAlignVariantCounts=as.data.frame(matrix(ncol=0,nrow=1))
f_KestrelAlignVariantCounts$AT_to_GC=sum(Fixed_wCpG$f_wCpG_AC) + sum(Fixed_wCpG$f_wCpG_AG) + sum(Fixed_wCpG$f_wCpG_TC) + sum(Fixed_wCpG$f_wCpG_TG)
f_KestrelAlignVariantCounts$GC_to_AT=sum(Fixed_wCpG$f_wCpG_GA) + sum(Fixed_wCpG$f_wCpG_GT) + sum(Fixed_wCpG$f_wCpG_CT) + sum(Fixed_wCpG$f_wCpG_CA)

f_KestrelAlignVariantCounts$AT_Total=sum(Fixed_wCpG$f_wCpG_AT) + sum(Fixed_wCpG$f_wCpG_AC) + sum(Fixed_wCpG$f_wCpG_AG) + sum(Fixed_wCpG$f_wCpG_TA) + sum(Fixed_wCpG$f_wCpG_TC) + sum(Fixed_wCpG$f_wCpG_TG)
f_KestrelAlignVariantCounts$GC_Total=sum(Fixed_wCpG$f_wCpG_GA) + sum(Fixed_wCpG$f_wCpG_GT) + sum(Fixed_wCpG$f_wCpG_GC) + sum(Fixed_wCpG$f_wCpG_CA) + sum(Fixed_wCpG$f_wCpG_CT) + sum(Fixed_wCpG$f_wCpG_CG)

f_KestrelAlignVariantCounts$AT_Equil=f_KestrelAlignVariantCounts$AT_Total - f_KestrelAlignVariantCounts$AT_to_GC + sum(AllrevIsochores100KB.ord$wCPG_CombinedAT)
f_KestrelAlignVariantCounts$GC_Equil=f_KestrelAlignVariantCounts$GC_Total - f_KestrelAlignVariantCounts$GC_to_AT + sum(AllrevIsochores100KB.ord$wCPG_CombinedGC)

#States chances of an equllibrium mutation are for GC and AT mutations--or that an AT site that mutates has a higher probability of staying AT than a GC site that mutates
Equillibruims=c(sum(f_KestrelAlignVariantCounts$AT_Equil), sum(f_KestrelAlignVariantCounts$GC_Equil))
Drifts=c(sum(f_KestrelAlignVariantCounts$AT_to_GC), sum(f_KestrelAlignVariantCounts$GC_to_AT))
Equillibrium_Drift.matrix=rbind(Equillibruims, Drifts)
colnames(Equillibrium_Drift.matrix)=c("AT", "GC")
chisq.test(Equillibrium_Drift.matrix)
OddsRatio=(Equillibrium_Drift.matrix[2,2]/Equillibrium_Drift.matrix[1,2])/(Equillibrium_Drift.matrix[2,1]/Equillibrium_Drift.matrix[1,1])
OddsRatio

f_wCpG_EquillibriumTable=data.frame(Genome=character(), ChiSquare=numeric(), Odds_Ratio=numeric(), P_Value=numeric())
for (GENOME in c("A", "B", "C", "D", "E", "F", "FpBSx77", "Kestrel")){
Equillibruims=c(sum(gFixed_wCpG[gFixed_wCpG$s==GENOME, c('f_wCpG_AT', 'f_wCpG_TA')], AllrevIsochores100KB.ord[AllrevIsochores100KB.ord$s==GENOME,'wCPG_CombinedAT']), sum(gFixed_wCpG[gFixed_wCpG$s==GENOME, c('f_wCpG_GC', 'f_wCpG_CG')], AllrevIsochores100KB.ord[AllrevIsochores100KB.ord$s==GENOME,'wCPG_CombinedGC']))
Drifts=c(sum(gFixed_wCpG[gFixed_wCpG$s==GENOME, c('f_wCpG_AC', 'f_wCpG_AG', 'f_wCpG_TC', 'f_wCpG_TG')]), sum(gFixed_wCpG[gFixed_wCpG$s==GENOME, c('f_wCpG_GA', 'f_wCpG_GT', 'f_wCpG_CT', 'f_wCpG_CA')]))
Equillibrium_Drift.matrix=rbind(Equillibruims, Drifts)
colnames(Equillibrium_Drift.matrix)=c("AT", "GC")
OddsRatio=(Equillibrium_Drift.matrix[2,2]/Equillibrium_Drift.matrix[1,2])/(Equillibrium_Drift.matrix[2,1]/Equillibrium_Drift.matrix[1,1])
CHISQ=chisq.test(Equillibrium_Drift.matrix)
ROW=data.frame(Genome=GENOME, ChiSquare=as.numeric(CHISQ$statistic), Odds_Ratio=as.numeric(OddsRatio), P_Value=as.numeric(CHISQ$p.value))
f_wCpG_EquillibriumTable=as.data.frame(rbind(f_wCpG_EquillibriumTable, ROW))
}

#Equillibria for heterozygous differences with CpG
h_KestrelAlignVariantCounts=as.data.frame(matrix(ncol=0,nrow=1))
h_KestrelAlignVariantCounts$AT_to_GC=sum(h_wCpG$h_wCpG_AC) + sum(h_wCpG$h_wCpG_AG) + sum(h_wCpG$h_wCpG_TC) + sum(h_wCpG$h_wCpG_TG)
h_KestrelAlignVariantCounts$GC_to_AT=sum(h_wCpG$h_wCpG_GA) + sum(h_wCpG$h_wCpG_GT) + sum(h_wCpG$h_wCpG_CT) + sum(h_wCpG$h_wCpG_CA)

h_KestrelAlignVariantCounts$AT_Total=sum(h_wCpG$h_wCpG_AT) + sum(h_wCpG$h_wCpG_AC) + sum(h_wCpG$h_wCpG_AG) + sum(h_wCpG$h_wCpG_TA) + sum(h_wCpG$h_wCpG_TC) + sum(h_wCpG$h_wCpG_TG)
h_KestrelAlignVariantCounts$GC_Total=sum(h_wCpG$h_wCpG_GA) + sum(h_wCpG$h_wCpG_GT) + sum(h_wCpG$h_wCpG_GC) + sum(h_wCpG$h_wCpG_CA) + sum(h_wCpG$h_wCpG_CT) + sum(h_wCpG$h_wCpG_CG)

h_KestrelAlignVariantCounts$AT_Equil=h_KestrelAlignVariantCounts$AT_Total - h_KestrelAlignVariantCounts$AT_to_GC + sum(AllrevIsochores100KB.ord$wCPG_CombinedAT)
h_KestrelAlignVariantCounts$GC_Equil=h_KestrelAlignVariantCounts$GC_Total - h_KestrelAlignVariantCounts$GC_to_AT + sum(AllrevIsochores100KB.ord$wCPG_CombinedGC)

#States chances of an equllibrium mutation are for GC and AT mutations--or that an AT site that mutates has a higher probability of staying AT than a GC site that mutates
Equillibruims=c(sum(h_KestrelAlignVariantCounts$AT_Equil), sum(h_KestrelAlignVariantCounts$GC_Equil))
Drifts=c(sum(h_KestrelAlignVariantCounts$AT_to_GC), sum(h_KestrelAlignVariantCounts$GC_to_AT))
Equillibrium_Drift.matrix=rbind(Equillibruims, Drifts)
colnames(Equillibrium_Drift.matrix)=c("AT", "GC")
chisq.test(Equillibrium_Drift.matrix)
OddsRatio=(Equillibrium_Drift.matrix[2,2]/Equillibrium_Drift.matrix[1,2])/(Equillibrium_Drift.matrix[2,1]/Equillibrium_Drift.matrix[1,1])
OddsRatio

h_wCpG_EquillibriumTable=data.frame(Genome=character(), ChiSquare=numeric(), Odds_Ratio=numeric(), P_Value=numeric())
for (GENOME in c("A", "B", "C", "D", "E", "F", "FpBSx77", "Kestrel")){
Equillibruims=c(sum(gh_wCpG[gh_wCpG$s==GENOME, c('h_wCpG_AT', 'h_wCpG_TA')], AllrevIsochores100KB.ord[AllrevIsochores100KB.ord$s==GENOME,'wCPG_CombinedAT']), sum(gh_wCpG[gh_wCpG$s==GENOME, c('h_wCpG_GC', 'h_wCpG_CG')], AllrevIsochores100KB.ord[AllrevIsochores100KB.ord$s==GENOME,'wCPG_CombinedGC']))
Drifts=c(sum(gh_wCpG[gh_wCpG$s==GENOME, c('h_wCpG_AC', 'h_wCpG_AG', 'h_wCpG_TC', 'h_wCpG_TG')]), sum(gh_wCpG[gh_wCpG$s==GENOME, c('h_wCpG_GA', 'h_wCpG_GT', 'h_wCpG_CT', 'h_wCpG_CA')]))
Equillibrium_Drift.matrix=rbind(Equillibruims, Drifts)
colnames(Equillibrium_Drift.matrix)=c("AT", "GC")
OddsRatio=(Equillibrium_Drift.matrix[2,2]/Equillibrium_Drift.matrix[1,2])/(Equillibrium_Drift.matrix[2,1]/Equillibrium_Drift.matrix[1,1])
CHISQ=chisq.test(Equillibrium_Drift.matrix)
ROW=data.frame(Genome=GENOME, ChiSquare=as.numeric(CHISQ$statistic), Odds_Ratio=as.numeric(OddsRatio), P_Value=as.numeric(CHISQ$p.value))
h_wCpG_EquillibriumTable=as.data.frame(rbind(h_wCpG_EquillibriumTable, ROW))
}

#Equillibria for fixed differences without CpG
f_KestrelAlignVariantCounts=as.data.frame(matrix(ncol=0,nrow=1))
f_KestrelAlignVariantCounts$AT_to_GC=sum(Fixed_nCpG$f_nCpG_AC) + sum(Fixed_nCpG$f_nCpG_AG) + sum(Fixed_nCpG$f_nCpG_TC) + sum(Fixed_nCpG$f_nCpG_TG)
f_KestrelAlignVariantCounts$GC_to_AT=sum(Fixed_nCpG$f_nCpG_GA) + sum(Fixed_nCpG$f_nCpG_GT) + sum(Fixed_nCpG$f_nCpG_CT) + sum(Fixed_nCpG$f_nCpG_CA)

f_KestrelAlignVariantCounts$AT_Total=sum(Fixed_nCpG$f_nCpG_AT) + sum(Fixed_nCpG$f_nCpG_AC) + sum(Fixed_nCpG$f_nCpG_AG) + sum(Fixed_nCpG$f_nCpG_TA) + sum(Fixed_nCpG$f_nCpG_TC) + sum(Fixed_nCpG$f_nCpG_TG)
f_KestrelAlignVariantCounts$GC_Total=sum(Fixed_nCpG$f_nCpG_GA) + sum(Fixed_nCpG$f_nCpG_GT) + sum(Fixed_nCpG$f_nCpG_GC) + sum(Fixed_nCpG$f_nCpG_CA) + sum(Fixed_nCpG$f_nCpG_CT) + sum(Fixed_nCpG$f_nCpG_CG)

f_KestrelAlignVariantCounts$AT_Equil=f_KestrelAlignVariantCounts$AT_Total - f_KestrelAlignVariantCounts$AT_to_GC + sum(na.omit(AllrevIsochores100KB.ord$nCpG_CombinedAT))
f_KestrelAlignVariantCounts$GC_Equil=f_KestrelAlignVariantCounts$GC_Total - f_KestrelAlignVariantCounts$GC_to_AT + sum(na.omit(AllrevIsochores100KB.ord$nCpG_CombinedGC))

#States chances of an equllibrium mutation are for GC and AT mutations--or that an AT site that mutates has a higher probability of staying AT than a GC site that mutates
Equillibruims=c(sum(f_KestrelAlignVariantCounts$AT_Equil), sum(f_KestrelAlignVariantCounts$GC_Equil))
Drifts=c(sum(f_KestrelAlignVariantCounts$AT_to_GC), sum(f_KestrelAlignVariantCounts$GC_to_AT))
Equillibrium_Drift.matrix=rbind(Equillibruims, Drifts)
colnames(Equillibrium_Drift.matrix)=c("AT", "GC")
chisq.test(Equillibrium_Drift.matrix)
OddsRatio=(Equillibrium_Drift.matrix[2,2]/Equillibrium_Drift.matrix[1,2])/(Equillibrium_Drift.matrix[2,1]/Equillibrium_Drift.matrix[1,1])
OddsRatio


f_nCpG_EquillibriumTable=data.frame(Genome=character(), ChiSquare=numeric(), Odds_Ratio=numeric(), P_Value=numeric())
for (GENOME in c("A", "B", "C", "D", "E", "F", "FpBSx77", "Kestrel")){
Equillibruims=c(sum(gFixed_nCpG[gFixed_nCpG$s==GENOME, c('f_nCpG_AT', 'f_nCpG_TA')], na.omit(AllrevIsochores100KB.ord[AllrevIsochores100KB.ord$s==GENOME,'nCpG_CombinedAT'])), sum(gFixed_nCpG[gFixed_nCpG$s==GENOME, c('f_nCpG_GC', 'f_nCpG_CG')], na.omit(AllrevIsochores100KB.ord[AllrevIsochores100KB.ord$s==GENOME,'nCpG_CombinedGC'])))
Drifts=c(sum(gFixed_nCpG[gFixed_nCpG$s==GENOME, c('f_nCpG_AC', 'f_nCpG_AG', 'f_nCpG_TC', 'f_nCpG_TG')]), sum(gFixed_nCpG[gFixed_nCpG$s==GENOME, c('f_nCpG_GA', 'f_nCpG_GT', 'f_nCpG_CT', 'f_nCpG_CA')]))
Equillibrium_Drift.matrix=rbind(Equillibruims, Drifts)
colnames(Equillibrium_Drift.matrix)=c("AT", "GC")
OddsRatio=(Equillibrium_Drift.matrix[2,2]/Equillibrium_Drift.matrix[1,2])/(Equillibrium_Drift.matrix[2,1]/Equillibrium_Drift.matrix[1,1])
CHISQ=chisq.test(Equillibrium_Drift.matrix)
ROW=data.frame(Genome=GENOME, ChiSquare=as.numeric(CHISQ$statistic), Odds_Ratio=as.numeric(OddsRatio), P_Value=as.numeric(CHISQ$p.value))
f_nCpG_EquillibriumTable=as.data.frame(rbind(f_nCpG_EquillibriumTable, ROW))
}

#Equillibria for heterozygous differences without CpG
h_KestrelAlignVariantCounts=as.data.frame(matrix(ncol=0,nrow=1))
h_KestrelAlignVariantCounts$AT_to_GC=sum(h_nCpG$h_nCpG_AC) + sum(h_nCpG$h_nCpG_AG) + sum(h_nCpG$h_nCpG_TC) + sum(h_nCpG$h_nCpG_TG)
h_KestrelAlignVariantCounts$GC_to_AT=sum(h_nCpG$h_nCpG_GA) + sum(h_nCpG$h_nCpG_GT) + sum(h_nCpG$h_nCpG_CT) + sum(h_nCpG$h_nCpG_CA)

h_KestrelAlignVariantCounts$AT_Total=sum(h_nCpG$h_nCpG_AT) + sum(h_nCpG$h_nCpG_AC) + sum(h_nCpG$h_nCpG_AG) + sum(h_nCpG$h_nCpG_TA) + sum(h_nCpG$h_nCpG_TC) + sum(h_nCpG$h_nCpG_TG)
h_KestrelAlignVariantCounts$GC_Total=sum(h_nCpG$h_nCpG_GA) + sum(h_nCpG$h_nCpG_GT) + sum(h_nCpG$h_nCpG_GC) + sum(h_nCpG$h_nCpG_CA) + sum(h_nCpG$h_nCpG_CT) + sum(h_nCpG$h_nCpG_CG)

h_KestrelAlignVariantCounts$AT_Equil=h_KestrelAlignVariantCounts$AT_Total - h_KestrelAlignVariantCounts$AT_to_GC + sum(na.omit(AllrevIsochores100KB.ord$nCpG_CombinedAT))
h_KestrelAlignVariantCounts$GC_Equil=h_KestrelAlignVariantCounts$GC_Total - h_KestrelAlignVariantCounts$GC_to_AT + sum(na.omit(AllrevIsochores100KB.ord$nCpG_CombinedGC))

#States chances of an equllibrium mutation are for GC and AT mutations--or that an AT site that mutates has a higher probability of staying AT than a GC site that mutates
Equillibruims=c(sum(h_KestrelAlignVariantCounts$AT_Equil), sum(h_KestrelAlignVariantCounts$GC_Equil))
Drifts=c(sum(h_KestrelAlignVariantCounts$AT_to_GC), sum(h_KestrelAlignVariantCounts$GC_to_AT))
Equillibrium_Drift.matrix=rbind(Equillibruims, Drifts)
colnames(Equillibrium_Drift.matrix)=c("AT", "GC")
chisq.test(Equillibrium_Drift.matrix)
OddsRatio=(Equillibrium_Drift.matrix[2,2]/Equillibrium_Drift.matrix[1,2])/(Equillibrium_Drift.matrix[2,1]/Equillibrium_Drift.matrix[1,1])
OddsRatio


h_nCpG_EquillibriumTable=data.frame(Genome=character(), ChiSquare=numeric(), Odds_Ratio=numeric(), P_Value=numeric())
for (GENOME in c("A", "B", "C", "D", "E", "F", "FpBSx77", "Kestrel")){
Equillibruims=c(sum(gh_nCpG[gh_nCpG$s==GENOME, c('h_nCpG_AT', 'h_nCpG_TA')], na.omit(AllrevIsochores100KB.ord[AllrevIsochores100KB.ord$s==GENOME,'nCpG_CombinedAT'])), sum(gh_nCpG[gh_nCpG$s==GENOME, c('h_nCpG_GC', 'h_nCpG_CG')], na.omit(AllrevIsochores100KB.ord[AllrevIsochores100KB.ord$s==GENOME,'nCpG_CombinedGC'])))
Drifts=c(sum(gh_nCpG[gh_nCpG$s==GENOME, c('h_nCpG_AC', 'h_nCpG_AG', 'h_nCpG_TC', 'h_nCpG_TG')]), sum(gh_nCpG[gh_nCpG$s==GENOME, c('h_nCpG_GA', 'h_nCpG_GT', 'h_nCpG_CT', 'h_nCpG_CA')]))
Equillibrium_Drift.matrix=rbind(Equillibruims, Drifts)
colnames(Equillibrium_Drift.matrix)=c("AT", "GC")
OddsRatio=(Equillibrium_Drift.matrix[2,2]/Equillibrium_Drift.matrix[1,2])/(Equillibrium_Drift.matrix[2,1]/Equillibrium_Drift.matrix[1,1])
CHISQ=chisq.test(Equillibrium_Drift.matrix)
ROW=data.frame(Genome=GENOME, ChiSquare=as.numeric(CHISQ$statistic), Odds_Ratio=as.numeric(OddsRatio), P_Value=as.numeric(CHISQ$p.value))
h_nCpG_EquillibriumTable=as.data.frame(rbind(h_nCpG_EquillibriumTable, ROW))
}


f_wCpG_EquillibriumTable$Type="w/CpG"
f_nCpG_EquillibriumTable$Type="n/CpG"
f_wCpG_EquillibriumTable$State="Fixed"
f_nCpG_EquillibriumTable$State="Fixed"
h_wCpG_EquillibriumTable$Type="w/CpG"
h_nCpG_EquillibriumTable$Type="n/CpG"
h_wCpG_EquillibriumTable$State="Het"
h_nCpG_EquillibriumTable$State="Het"
f_wCpG_EquillibriumTable$Group="Fixed w/CpG"
f_nCpG_EquillibriumTable$Group="Fixed n/CpG"
h_wCpG_EquillibriumTable$Group="Het. w/CpG"
h_nCpG_EquillibriumTable$Group="Het. n/CpG"
h_nCpG_EquillibriumTable$Group=factor(h_nCpG_EquillibriumTable$Group, levels=c("Het. w/CpG","Fixed w/CpG","Het. n/CpG","Fixed n/CpG"))


AllEqullibriumTable=rbind(f_wCpG_EquillibriumTable,f_nCpG_EquillibriumTable,h_wCpG_EquillibriumTable,h_nCpG_EquillibriumTable)
AllEqullibriumTable=AllEqullibriumTable[AllEqullibriumTable$Odds_Ratio!="NaN",]

summary(aov(AllEqullibriumTable$Odds_Ratio~AllEqullibriumTable$Type*AllEqullibriumTable$State))
pairwise.t.test(AllEqullibriumTable$Odds_Ratio, AllEqullibriumTable$Group)
aggregate(AllEqullibriumTable$Odds_Ratio, list(AllEqullibriumTable$Group), FUN=mean)

Graph_EquillibriumTable=AllEqullibriumTable
Graph_EquillibriumTable$Genome=as.character(Graph_EquillibriumTable$Genome)
Graph_EquillibriumTable$Genome[Graph_EquillibriumTable$Genome=="A"]="Lanner"
Graph_EquillibriumTable$Genome[Graph_EquillibriumTable$Genome=="B"]="Barbary"
Graph_EquillibriumTable$Genome[Graph_EquillibriumTable$Genome=="C"]="Peregrine"
Graph_EquillibriumTable$Genome[Graph_EquillibriumTable$Genome=="D"]="Saker"
Graph_EquillibriumTable$Genome[Graph_EquillibriumTable$Genome=="E"]="Gyr-1"
Graph_EquillibriumTable$Genome[Graph_EquillibriumTable$Genome=="F"]="Gyr-2"
Graph_EquillibriumTable$Genome[Graph_EquillibriumTable$Genome=="FpBSx77"]="Bl. Shaheen"
Graph_EquillibriumTable$Genome[Graph_EquillibriumTable$Genome=="Kestrel"]="C. Kestrel"

Graph_EquillibriumTable$Equillibrium_GC=(1-(Graph_EquillibriumTable$Odds_Ratio/(Graph_EquillibriumTable$Odds_Ratio+1)))*100

Export_EquillibriumTable=Graph_EquillibriumTable[Graph_EquillibriumTable$Type=="w/CpG",]
Export_EquillibriumTable=Export_EquillibriumTable[Export_EquillibriumTable$State=="Fixed",]
Het_EquillibriumTable=Graph_EquillibriumTable[Graph_EquillibriumTable$Type=="w/CpG",]
Het_EquillibriumTable=Het_EquillibriumTable[Het_EquillibriumTable$State=="Het",]
Het_EquillibriumTable=Het_EquillibriumTable[Het_EquillibriumTable$Genome!="C. Kestrel",]
Het_EquillibriumTable.ord=Het_EquillibriumTable[order(Het_EquillibriumTable$Genome),]
setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis")
SummarySupernova=read.delim("Summary_SuperNovaAssembliesRef8.txt", header=TRUE, row.names=1, sep="\t")
SummarySupernova$Genome=row.names(SummarySupernova)
SummarySupernova=SummarySupernova[SummarySupernova$Genome!="",]

SummarySupernova$Genome=as.character(SummarySupernova$Genome)
SummarySupernova$Genome[SummarySupernova$Genome=="A"]="Lanner"
SummarySupernova$Genome[SummarySupernova$Genome=="B"]="Barbary"
SummarySupernova$Genome[SummarySupernova$Genome=="C"]="Peregrine"
SummarySupernova$Genome[SummarySupernova$Genome=="D"]="Saker"
SummarySupernova$Genome[SummarySupernova$Genome=="E"]="Gyr-1"
SummarySupernova$Genome[SummarySupernova$Genome=="F"]="Gyr-2"
SummarySupernova$Genome[SummarySupernova$Genome=="FpBSx77"]="Bl. Shaheen"
SummarySupernova$Genome[SummarySupernova$Genome=="Kestrel"]="C. Kestrel"
SummarySupernova.purge=SummarySupernova[SummarySupernova$Genome!="C. Kestrel",]

SummarySupernova.ord=SummarySupernova.purge[order(SummarySupernova.purge$Genome),]
Export_EquillibriumTable.ord=Export_EquillibriumTable[order(Export_EquillibriumTable$Genome),]
Export_EquillibriumTable.ord$HetEquillibrium=Het_EquillibriumTable.ord$Equillibrium
Export_EquillibriumTable.ord$pGC=SummarySupernova.ord$Assembly_GC_Content

setwd("~/NYU_AD/FalconProject/Manuscripts/FalconRef8/UltimateTables/SupplementaryTables/Intermediate")
write.table(Export_EquillibriumTable.ord, file="wCpG_SuppTable_Equillibrium.txt", sep="\t", row.names=FALSE)

Export_EquillibriumTable=Graph_EquillibriumTable[Graph_EquillibriumTable$Type=="n/CpG",]
Export_EquillibriumTable=Export_EquillibriumTable[Export_EquillibriumTable$State=="Fixed",]
Het_EquillibriumTable=Graph_EquillibriumTable[Graph_EquillibriumTable$Type=="n/CpG",]
Het_EquillibriumTable=Het_EquillibriumTable[Het_EquillibriumTable$State=="Het",]
Het_EquillibriumTable=Het_EquillibriumTable[Het_EquillibriumTable$Genome!="C. Kestrel",]
Het_EquillibriumTable.ord=Het_EquillibriumTable[order(Het_EquillibriumTable$Genome),]
Export_EquillibriumTable.ord=Export_EquillibriumTable[order(Export_EquillibriumTable$Genome),]
Export_EquillibriumTable.ord$HetEquillibrium=Het_EquillibriumTable.ord$Equillibrium
Export_EquillibriumTable.ord$pGC=SummarySupernova.ord$Assembly_GC_Content
write.table(Export_EquillibriumTable.ord, file="nCpG_SuppTable_Equillibrium.txt", sep="\t", row.names=FALSE)


library(ggplot2)
ggplot(Graph_EquillibriumTable, aes(x=Genome, y=Odds_Ratio, group=Group)) + geom_point(aes(color=Group, shape=Group), size=5) + theme(text= element_text(size=12)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, color = "black", face="bold", size=12), legend.text=element_text(size=10), legend.title = element_blank()) + labs(x="",y="Odds Ratio G/C to A/T") + geom_hline(yintercept=1, linetype="dashed", color = "red", size=2) + ylim(0,1.8) + scale_colour_manual(values=c('Fixed w/CpG'="blue", 'Het. w/CpG'="blue", 'Het. n/CpG'="black", 'Fixed n/CpG'="black")) + scale_linetype_manual(values=c('Fixed w/CpG'="solid", 'Het. w/CpG'="dotted", 'Het. n/CpG'="dotted", 'Fixed n/CpG'="solid")) + scale_shape_manual(values=c('Fixed w/CpG'=15, 'Het. w/CpG'=16, 'Het. n/CpG'=17, 'Fixed n/CpG'=18)) + theme(legend.position = c(0.85, 0.25))
ggsave("suppFig6B_AT-GC_Equillibrium.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Intermediate", scale = 1, width = 3.6, height = 5, units = c("in"), dpi =1200, limitsize = TRUE)

setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis/Synthesis/Data/ChromAnnotated")
Isochores100KB.unord=read.delim("AllultChromAnnotated100KBWindowsSummary.txt", header=TRUE, sep="\t")
Isochores100KB.nofix=Isochores100KB.unord[order(Isochores100KB.unord$Genome,Isochores100KB.unord$Window),]
Isochores100KB=cbind(Isochores100KB.nofix,Fixed_wCpG,Fixed_nCpG,h_wCpG,h_nCpG)

rm(Isochores100KB.unord)
rm(Isochores100KB.nofix)

Isochores100KB$Genome=as.character(Isochores100KB$Genome)
Isochores100KB$Genome[Isochores100KB$Genome=="A"]="Lanner"
Isochores100KB$Genome[Isochores100KB$Genome=="B"]="Barbary"
Isochores100KB$Genome[Isochores100KB$Genome=="C"]="Peregrine"
Isochores100KB$Genome[Isochores100KB$Genome=="D"]="Saker"
Isochores100KB$Genome[Isochores100KB$Genome=="E"]="Gyr-1"
Isochores100KB$Genome[Isochores100KB$Genome=="F"]="Gyr-2"
Isochores100KB$Genome[Isochores100KB$Genome=="FpBSx77"]="Bl. Shaheen"
Isochores100KB$Genome[Isochores100KB$Genome=="Kestrel"]="C. Kestrel"

Isochores100KB$Family=rep("NA", length(Isochores100KB$pGC))
Isochores100KB$Family[Isochores100KB$pGC >=0 & Isochores100KB$pGC<0.37]="L1"
Isochores100KB$Family[Isochores100KB$pGC>=0.37 & Isochores100KB$pGC<0.41]="L2"
Isochores100KB$Family[Isochores100KB$pGC>=0.41 & Isochores100KB$pGC<0.46]="H1"
Isochores100KB$Family[Isochores100KB$pGC>=0.46 & Isochores100KB$pGC<0.53]="H2"
Isochores100KB$Family[Isochores100KB$pGC>=0.53 & Isochores100KB$pGC<=1]="H3"
Isochores100KB$Family=factor(Isochores100KB$Family, levels=c("L1","L2","H1","H2","H3"))
library(ggplot2)
mycolors=c("#051094","#48AAAD","#EFFD5F","#FC7F03","#DC143C")
Isochores100KB$Genome=factor(Isochores100KB$Genome, levels = c("Barbary", "Bl. Shaheen", "Peregrine", "Lanner", "Gyr-1", "Gyr-2", "Saker", "C. Kestrel"))
ggplot(Isochores100KB, aes(x=pGC, fill=Family, color=Family)) + geom_histogram(binwidth = 0.005, breaks=seq(.34, .70, 0.005), aes(weight=tBASES/1000000), size=0.25) + facet_wrap(vars(Genome), nrow=2, ncol=4) + labs(x="Percent GC", y="Megabases") + scale_color_manual(values=rep("black", 5)) + scale_fill_manual(values=mycolors) + theme(strip.text=element_text(face="bold"))
ggsave("suppFig3_FalcoIsochores100KB.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Supplementary", scale = 1, width = 7, height = 5, units = c("in"), dpi =1200, limitsize = TRUE)
RepIsochores100KB=Isochores100KB[Isochores100KB$Genome=="Lanner",c(1:12,130)]

Isochores1MB=read.delim("AllultChromAnnotated1MBWindowsSummary.txt", header=TRUE, sep="\t")
Isochores1MB$Genome=as.character(Isochores1MB$Genome)
Isochores1MB$Genome[Isochores1MB$Genome=="A"]="Lanner"
Isochores1MB$Genome[Isochores1MB$Genome=="B"]="Barbary"
Isochores1MB$Genome[Isochores1MB$Genome=="C"]="Peregrine"
Isochores1MB$Genome[Isochores1MB$Genome=="D"]="Saker"
Isochores1MB$Genome[Isochores1MB$Genome=="E"]="Gyr-1"
Isochores1MB$Genome[Isochores1MB$Genome=="F"]="Gyr-2"
Isochores1MB$Genome[Isochores1MB$Genome=="FpBSx77"]="Bl. Shaheen"
Isochores1MB$Genome[Isochores1MB$Genome=="Kestrel"]="C. Kestrel"

Isochores1MB$Family=rep("NA", length(Isochores1MB$pGC))
Isochores1MB$Family[Isochores1MB$pGC >=0 & Isochores1MB$pGC<0.37]="L1"
Isochores1MB$Family[Isochores1MB$pGC>=0.37 & Isochores1MB$pGC<0.41]="L2"
Isochores1MB$Family[Isochores1MB$pGC>=0.41 & Isochores1MB$pGC<0.46]="H1"
Isochores1MB$Family[Isochores1MB$pGC>=0.46 & Isochores1MB$pGC<0.53]="H2"
Isochores1MB$Family[Isochores1MB$pGC>=0.53 & Isochores1MB$pGC<=1]="H3"
Isochores1MB$Family=factor(Isochores1MB$Family, levels=c("L1","L2","H1","H2","H3"))
library(ggplot2)
mycolors=c("#051094","#48AAAD","#EFFD5F","#FC7F03","#DC143C")
Isochores1MB$Genome=factor(Isochores1MB$Genome, levels = c("Barbary", "Bl. Shaheen", "Peregrine", "Lanner", "Gyr-1", "Gyr-2", "Saker", "C. Kestrel"))
ggplot(Isochores1MB, aes(x=pGC, fill=Family, color=Family)) + geom_histogram(binwidth = 0.005, breaks=seq(.34, .70, 0.005), aes(weight=tBASES/1000000), size=0.25) + facet_wrap(vars(Genome), nrow=2, ncol=4) + labs(x="Percent GC", y="Megabases") + scale_color_manual(values=rep("black", 5)) + scale_fill_manual(values=mycolors) + theme(strip.text=element_text(face="bold"))
ggsave("suppFig4_FalcoIsochores1MB.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Supplementary/", scale = 1, width = 7, height = 5, units = c("in"), dpi =1200, limitsize = TRUE)

#R
setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis/Synthesis/Data/OtherGenomes_Full")
IsochoresOther=read.delim("AllOtherGenomeFull100KBWindowsSummary.txt", header=TRUE, sep="\t")

IsochoresOther$Family=rep("NA", length(IsochoresOther$pGC))
IsochoresOther$Family[IsochoresOther$pGC >=0 & IsochoresOther$pGC<0.37]="L1"
IsochoresOther$Family[IsochoresOther$pGC>=0.37 & IsochoresOther$pGC<0.41]="L2"
IsochoresOther$Family[IsochoresOther$pGC>=0.41 & IsochoresOther$pGC<0.46]="H1"
IsochoresOther$Family[IsochoresOther$pGC>=0.46 & IsochoresOther$pGC<0.53]="H2"
IsochoresOther$Family[IsochoresOther$pGC>=0.53 & IsochoresOther$pGC<=1]="H3"
IsochoresOther$Family=factor(IsochoresOther$Family, levels=c("L1","L2","H1","H2","H3"))
library(ggplot2)
mycolors=c("#051094","#48AAAD","#EFFD5F","#FC7F03","#DC143C")
IsochoresOther$Genome=as.character(IsochoresOther$Genome)
IsochoresOther$Genome[IsochoresOther$Genome=="Gallus_gallus"]="Gallus"
IsochoresOther$Genome=factor(IsochoresOther$Genome, levels = c("Human", "Lacerta", "Crocodile", "Turtle", "Gallus", "Seriema", "Kakapo", "Thrush"))
ggplot(IsochoresOther, aes(x=pGC, fill=Family, color=Family)) + geom_histogram(binwidth = 0.005, breaks=seq(.34, .70, 0.005), aes(weight=tBASES/1000000), size=0.25) + facet_wrap(vars(Genome), nrow=2, ncol=4) + labs(x="Percent GC", y="Megabases") + scale_color_manual(values=rep("black", 5)) + scale_fill_manual(values=mycolors) + theme(strip.text=element_text(face="bold"))
ggsave("Unused_OtherGenomes100KB.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Supplementary/Unused", scale = 1, width = 7, height = 5, units = c("in"), dpi =1200, limitsize = TRUE)

RepIsochores100KB=as.data.frame(rbind(RepIsochores100KB, IsochoresOther[IsochoresOther$Genome=="Human",], IsochoresOther[IsochoresOther$Genome=="Lacerta",], IsochoresOther[IsochoresOther$Genome=="Turtle",], IsochoresOther[IsochoresOther$Genome=="Crocodile",], IsochoresOther[IsochoresOther$Genome=="Thrush",], IsochoresOther[IsochoresOther$Genome=="Seriema",], IsochoresOther[IsochoresOther$Genome=="Kakapo",]))
RepIsochores100KB$Genome=factor(RepIsochores100KB$Genome, levels = c("Lanner", "Kakapo", "Thrush", "Seriema", "Crocodile", "Turtle", "Lacerta", "Human"))
ggplot(RepIsochores100KB, aes(x=pGC, fill=Family, color=Family)) + geom_histogram(binwidth = 0.005, breaks=seq(.34, .70, 0.005), aes(weight=tBASES/1000000), size=0.25) + facet_wrap(vars(Genome), nrow=2, ncol=4) + labs(x="Percent GC", y="Megabases") + scale_color_manual(values=rep("black", 5)) + scale_fill_manual(values=mycolors) + theme(strip.text=element_text(face="bold"))
ggsave("Fig3_RepIsochores100KB.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/", scale = 1, width = 7, height = 5, units = c("in"), dpi =1200, limitsize = TRUE)

IsochoresOther=read.delim("AllOtherGenomeFull1MBWindowsSummary.txt", header=TRUE, sep="\t")

IsochoresOther$Family=rep("NA", length(IsochoresOther$pGC))
IsochoresOther$Family[IsochoresOther$pGC >=0 & IsochoresOther$pGC<0.37]="L1"
IsochoresOther$Family[IsochoresOther$pGC>=0.37 & IsochoresOther$pGC<0.41]="L2"
IsochoresOther$Family[IsochoresOther$pGC>=0.41 & IsochoresOther$pGC<0.46]="H1"
IsochoresOther$Family[IsochoresOther$pGC>=0.46 & IsochoresOther$pGC<0.53]="H2"
IsochoresOther$Family[IsochoresOther$pGC>=0.53 & IsochoresOther$pGC<=1]="H3"
IsochoresOther$Family=factor(IsochoresOther$Family, levels=c("L1","L2","H1","H2","H3"))
library(ggplot2)
mycolors=c("#051094","#48AAAD","#EFFD5F","#FC7F03","#DC143C")
IsochoresOther$Genome=as.character(IsochoresOther$Genome)
IsochoresOther$Genome[IsochoresOther$Genome=="Gallus_gallus"]="Gallus"
IsochoresOther$Genome=factor(IsochoresOther$Genome, levels = c("Human", "Lacerta", "Crocodile", "Turtle", "Gallus", "Seriema", "Kakapo", "Thrush"))
ggplot(IsochoresOther, aes(x=pGC, fill=Family, color=Family)) + geom_histogram(binwidth = 0.005, breaks=seq(.34, .70, 0.005), aes(weight=tBASES/1000000), size=0.25) + facet_wrap(vars(Genome), nrow=2, ncol=4) + labs(x="Percent GC", y="Megabases") + scale_color_manual(values=rep("black", 5)) + scale_fill_manual(values=mycolors) + theme(strip.text=element_text(face="bold"))
ggsave("Unused_OtherGenomes1MB.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Supplementary/Unused", scale = 1, width = 7, height = 5, units = c("in"), dpi =1200, limitsize = TRUE)

Isochores100KB$f_AT_wCpG_Equillibirum=Isochores100KB$f_wCpG_AT+Isochores100KB$f_wCpG_TA+Isochores100KB$wCPG_CombinedAT
Isochores100KB$f_AT_nCpG_Equillibirum=Isochores100KB$f_nCpG_AT+Isochores100KB$f_nCpG_TA+Isochores100KB$nCpG_CombinedAT

Isochores100KB$f_GC_wCpG_Equillibirum=Isochores100KB$f_wCpG_GC+Isochores100KB$f_wCpG_CG+Isochores100KB$wCPG_CombinedGC
Isochores100KB$f_GC_nCpG_Equillibirum=Isochores100KB$f_nCpG_GC+Isochores100KB$f_nCpG_CG+Isochores100KB$nCpG_CombinedGC

Isochores100KB$f_ATtoGC_wCpG_Drifts=Isochores100KB$f_wCpG_AG+Isochores100KB$f_wCpG_AC+Isochores100KB$f_wCpG_TG+Isochores100KB$f_wCpG_TC
Isochores100KB$f_ATtoGC_nCpG_Drifts=Isochores100KB$f_nCpG_AG+Isochores100KB$f_nCpG_AC+Isochores100KB$f_nCpG_TG+Isochores100KB$f_nCpG_TC

Isochores100KB$f_GCtoAT_wCpG_Drifts=Isochores100KB$f_wCpG_GA+Isochores100KB$f_wCpG_CA+Isochores100KB$f_wCpG_GT+Isochores100KB$f_wCpG_CT
Isochores100KB$f_GCtoAT_nCpG_Drifts=Isochores100KB$f_nCpG_GA+Isochores100KB$f_nCpG_CA+Isochores100KB$f_nCpG_GT+Isochores100KB$f_nCpG_CT

Isochores100KB$f_wCpG_OddsRatioGCtoAT=(Isochores100KB$f_GCtoAT_wCpG_Drifts+1/Isochores100KB$f_GC_wCpG_Equillibirum)/(Isochores100KB$f_ATtoGC_wCpG_Drifts+1/Isochores100KB$f_AT_wCpG_Equillibirum)
Isochores100KB$f_nCpG_OddsRatioGCtoAT=(Isochores100KB$f_GCtoAT_nCpG_Drifts+1/Isochores100KB$f_GC_nCpG_Equillibirum)/(Isochores100KB$f_ATtoGC_nCpG_Drifts+1/Isochores100KB$f_AT_nCpG_Equillibirum)

Isochores100KB$h_AT_wCpG_Equillibirum=Isochores100KB$h_wCpG_AT+Isochores100KB$h_wCpG_TA+Isochores100KB$wCPG_CombinedAT
Isochores100KB$h_AT_nCpG_Equillibirum=Isochores100KB$h_nCpG_AT+Isochores100KB$h_nCpG_TA+Isochores100KB$nCpG_CombinedAT

Isochores100KB$h_GC_wCpG_Equillibirum=Isochores100KB$h_wCpG_GC+Isochores100KB$h_wCpG_CG+Isochores100KB$wCPG_CombinedGC
Isochores100KB$h_GC_nCpG_Equillibirum=Isochores100KB$h_nCpG_GC+Isochores100KB$h_nCpG_CG+Isochores100KB$nCpG_CombinedGC

Isochores100KB$h_ATtoGC_wCpG_Drifts=Isochores100KB$h_wCpG_AG+Isochores100KB$h_wCpG_AC+Isochores100KB$h_wCpG_TG+Isochores100KB$h_wCpG_TC
Isochores100KB$h_ATtoGC_nCpG_Drifts=Isochores100KB$h_nCpG_AG+Isochores100KB$h_nCpG_AC+Isochores100KB$h_nCpG_TG+Isochores100KB$h_nCpG_TC

Isochores100KB$h_GCtoAT_wCpG_Drifts=Isochores100KB$h_wCpG_GA+Isochores100KB$h_wCpG_CA+Isochores100KB$h_wCpG_GT+Isochores100KB$h_wCpG_CT
Isochores100KB$h_GCtoAT_nCpG_Drifts=Isochores100KB$h_nCpG_GA+Isochores100KB$h_nCpG_CA+Isochores100KB$h_nCpG_GT+Isochores100KB$h_nCpG_CT

Isochores100KB$h_wCpG_OddsRatioGCtoAT=(Isochores100KB$h_GCtoAT_wCpG_Drifts+1/Isochores100KB$h_GC_wCpG_Equillibirum)/(Isochores100KB$h_ATtoGC_wCpG_Drifts+1/Isochores100KB$h_AT_wCpG_Equillibirum)
Isochores100KB$h_nCpG_OddsRatioGCtoAT=(Isochores100KB$h_GCtoAT_nCpG_Drifts+1/Isochores100KB$h_GC_nCpG_Equillibirum)/(Isochores100KB$h_ATtoGC_nCpG_Drifts+1/Isochores100KB$h_AT_nCpG_Equillibirum)

Isochores100KB$t_AT_wCpG_Equillibirum=Isochores100KB$f_wCpG_AT+Isochores100KB$f_wCpG_TA+Isochores100KB$h_wCpG_AT+Isochores100KB$h_wCpG_TA+Isochores100KB$wCPG_CombinedAT
Isochores100KB$t_AT_nCpG_Equillibirum=Isochores100KB$f_nCpG_AT+Isochores100KB$f_nCpG_TA+Isochores100KB$h_nCpG_AT+Isochores100KB$h_nCpG_TA+Isochores100KB$nCpG_CombinedAT

Isochores100KB$t_GC_wCpG_Equillibirum=Isochores100KB$f_wCpG_GC+Isochores100KB$f_wCpG_CG+Isochores100KB$h_wCpG_GC+Isochores100KB$h_wCpG_CG+Isochores100KB$wCPG_CombinedGC
Isochores100KB$t_GC_nCpG_Equillibirum=Isochores100KB$f_nCpG_GC+Isochores100KB$f_nCpG_CG+Isochores100KB$h_nCpG_GC+Isochores100KB$h_nCpG_CG+Isochores100KB$nCpG_CombinedGC

Isochores100KB$t_ATtoGC_wCpG_Drifts=Isochores100KB$f_wCpG_AG+Isochores100KB$f_wCpG_AC+Isochores100KB$f_wCpG_TG+Isochores100KB$f_wCpG_TC+Isochores100KB$h_wCpG_AG+Isochores100KB$h_wCpG_AC+Isochores100KB$h_wCpG_TG+Isochores100KB$h_wCpG_TC
Isochores100KB$t_ATtoGC_nCpG_Drifts=Isochores100KB$f_nCpG_AG+Isochores100KB$f_nCpG_AC+Isochores100KB$f_nCpG_TG+Isochores100KB$f_nCpG_TC+Isochores100KB$h_nCpG_AG+Isochores100KB$h_nCpG_AC+Isochores100KB$h_nCpG_TG+Isochores100KB$h_nCpG_TC

Isochores100KB$t_GCtoAT_wCpG_Drifts=Isochores100KB$f_wCpG_GA+Isochores100KB$f_wCpG_CA+Isochores100KB$f_wCpG_GT+Isochores100KB$f_wCpG_CT+Isochores100KB$h_wCpG_GA+Isochores100KB$h_wCpG_CA+Isochores100KB$h_wCpG_GT+Isochores100KB$h_wCpG_CT    
Isochores100KB$t_GCtoAT_nCpG_Drifts=Isochores100KB$f_nCpG_GA+Isochores100KB$f_nCpG_CA+Isochores100KB$f_nCpG_GT+Isochores100KB$f_nCpG_CT+Isochores100KB$h_nCpG_GA+Isochores100KB$h_nCpG_CA+Isochores100KB$h_nCpG_GT+Isochores100KB$h_nCpG_CT

Isochores100KB$t_wCpG_OddsRatioGCtoAT=(Isochores100KB$t_GCtoAT_wCpG_Drifts/Isochores100KB$t_GC_wCpG_Equillibirum+1)/(Isochores100KB$t_ATtoGC_wCpG_Drifts+1/Isochores100KB$t_AT_wCpG_Equillibirum)
Isochores100KB$t_nCpG_OddsRatioGCtoAT=(Isochores100KB$t_GCtoAT_nCpG_Drifts/Isochores100KB$t_GC_nCpG_Equillibirum+1)/(Isochores100KB$t_ATtoGC_nCpG_Drifts+1/Isochores100KB$t_AT_nCpG_Equillibirum)

Isochores100KB_wCpG=Isochores100KB[is.finite(Isochores100KB$f_wCpG_OddsRatioGCtoAT),]
Isochores100KB_nCpG=Isochores100KB[is.finite(Isochores100KB$f_nCpG_OddsRatioGCtoAT),]
Isochores100KB_wCpG=Isochores100KB_wCpG[Isochores100KB_wCpG$f_wCpG_OddsRatioGCtoAT>0,]
Isochores100KB_nCpG=Isochores100KB_nCpG[Isochores100KB_nCpG$f_nCpG_OddsRatioGCtoAT>0,]

Isochores100KB_wCpG.noamb=Isochores100KB_wCpG[Isochores100KB_wCpG$FalcoCONSENSUS_CLASS_MATCH!="Ambiguous",]
Isochores100KB_wCpG.noamb=Isochores100KB_wCpG.noamb[Isochores100KB_wCpG.noamb$CONSENSUS_CLASS_MATCH!="Ambiguous",]

Isochores100KB_wCpG.noamb$ANCESTRAL_STATE=rep("NA", length(Isochores100KB_wCpG.noamb$CONSENSUS_CLASS_MATCH))
Isochores100KB_wCpG.noamb$ANCESTRAL_STATE[Isochores100KB_wCpG.noamb$CONSENSUS_CLASS_MATCH!="Microchromosome"]="Large"
Isochores100KB_wCpG.noamb$ANCESTRAL_STATE[Isochores100KB_wCpG.noamb$CONSENSUS_CLASS_MATCH=="Microchromosome"]="Microchromosome"

Isochores100KB_wCpG.noamb$FALCO_STATE=rep("NA", length(Isochores100KB_wCpG.noamb$FalcoCONSENSUS_CLASS_MATCH))
Isochores100KB_wCpG.noamb$FALCO_STATE[Isochores100KB_wCpG.noamb$FalcoCONSENSUS_CLASS_MATCH!="Microchromosome"]="Large"
Isochores100KB_wCpG.noamb$FALCO_STATE[Isochores100KB_wCpG.noamb$FalcoCONSENSUS_CLASS_MATCH=="Microchromosome"]="Microchromosome"

Isochores100KB_wCpG.noamb$ChromosomeState=paste(Isochores100KB_wCpG.noamb$ANCESTRAL_STATE, Isochores100KB_wCpG.noamb$FALCO_STATE, sep=" | ")
Isochores100KB_wCpG.noamb$ChromosomeState=gsub("Microchromosome","Micro",Isochores100KB_wCpG.noamb$ChromosomeState)
Isochores100KB_wCpG.noamb$ChromosomeState=factor(Isochores100KB_wCpG.noamb$ChromosomeState, c("Large | Micro", "Large | Large", "Micro | Large", "Micro | Micro"))

Isochores100KB_nCpG.noamb=Isochores100KB_nCpG[Isochores100KB_nCpG$FalcoCONSENSUS_CLASS_MATCH!="Ambiguous",]
Isochores100KB_nCpG.noamb=Isochores100KB_nCpG.noamb[Isochores100KB_nCpG.noamb$CONSENSUS_CLASS_MATCH!="Ambiguous",]

Isochores100KB_nCpG.noamb$ANCESTRAL_STATE=rep("NA", length(Isochores100KB_nCpG.noamb$CONSENSUS_CLASS_MATCH))
Isochores100KB_nCpG.noamb$ANCESTRAL_STATE[Isochores100KB_nCpG.noamb$CONSENSUS_CLASS_MATCH!="Microchromosome"]="Large"
Isochores100KB_nCpG.noamb$ANCESTRAL_STATE[Isochores100KB_nCpG.noamb$CONSENSUS_CLASS_MATCH=="Microchromosome"]="Microchromosome"

Isochores100KB_nCpG.noamb$FALCO_STATE=rep("NA", length(Isochores100KB_nCpG.noamb$FalcoCONSENSUS_CLASS_MATCH))
Isochores100KB_nCpG.noamb$FALCO_STATE[Isochores100KB_nCpG.noamb$FalcoCONSENSUS_CLASS_MATCH!="Microchromosome"]="Large"
Isochores100KB_nCpG.noamb$FALCO_STATE[Isochores100KB_nCpG.noamb$FalcoCONSENSUS_CLASS_MATCH=="Microchromosome"]="Microchromosome"

Isochores100KB_nCpG.noamb$ChromosomeState=paste(Isochores100KB_nCpG.noamb$ANCESTRAL_STATE, Isochores100KB_nCpG.noamb$FALCO_STATE, sep=" | ")
Isochores100KB_nCpG.noamb$ChromosomeState=gsub("Microchromosome","Micro",Isochores100KB_nCpG.noamb$ChromosomeState)
Isochores100KB_nCpG.noamb$ChromosomeState=factor(Isochores100KB_nCpG.noamb$ChromosomeState, c("Large | Micro", "Large | Large", "Micro | Large", "Micro | Micro"))

library(multcompView)
library(nortest)
setwd("~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Supplementary/")

#Chromosome State Differences Test for Mutation Odds Ration w/CpG
f_wCPG.anova.nocontgc=aov(log2(f_wCpG_OddsRatioGCtoAT+1)~Genome+ANCESTRAL_STATE:FALCO_STATE, data=Isochores100KB_wCpG.noamb)
f_wCpGnocontgc_TUKEY=TukeyHSD(x=f_wCPG.anova.nocontgc, 'ANCESTRAL_STATE:FALCO_STATE', conf.level=0.95)
summary(f_wCPG.anova.nocontgc)
f_wCpGnocontgc_TUKEY
plot(f_wCpGnocontgc_TUKEY, las=1 , col="brown", cex=0.5)
png("suppFig11f_wCpGnocontgc_TUKEY.png", width=14, height=7, units="in", res=1200)
par(mar=c(5,26,4,1)+.1)
plot(f_wCpGnocontgc_TUKEY, las=1 , col="brown", cex=0.5)
dev.off()

f_wCPG.anova.gc=aov(log2(f_wCpG_OddsRatioGCtoAT+1)~pGC, data=Isochores100KB_wCpG.noamb)
f_wCPG.anova.gcres=residuals(f_wCPG.anova.gc)
ad.test(f_wCPG.anova.gcres)
f_wCPG.anova=aov(f_wCPG.anova.gcres~Genome+ANCESTRAL_STATE:FALCO_STATE, data=Isochores100KB_wCpG.noamb)
summary(f_wCPG.anova.gc)
summary.lm(f_wCPG.anova.gc)
summary(f_wCPG.anova)

f_nCPG.anova.gc=aov(log2(f_nCpG_OddsRatioGCtoAT+1)~pGC, data=Isochores100KB_nCpG.noamb)
f_nCPG.anova.gcres=residuals(f_nCPG.anova.gc)
ad.test(f_nCPG.anova.gcres)
f_nCPG.anova=aov(f_nCPG.anova.gcres~Genome+ANCESTRAL_STATE:FALCO_STATE, data=Isochores100KB_nCpG.noamb)
summary(f_nCPG.anova.gc)
summary.lm(f_nCPG.anova.gc)
summary(f_nCPG.anova)

f_wCPG.anova.gc=aov(log2(f_wCpG_OddsRatioGCtoAT+1)~pGC, data=Isochores100KB_wCpG.noamb)
f_wCPG.anova.gcres=residuals(f_wCPG.anova.gc)
ad.test(f_wCPG.anova.gcres)
f_wCPG.anova=aov(f_wCPG.anova.gcres~Genome+ANCESTRAL_STATE:FALCO_STATE, data=Isochores100KB_wCpG.noamb)
summary(f_wCPG.anova.gc)
summary(f_wCPG.anova)

f_wCpG_BoxPlot=ggplot(Isochores100KB_wCpG.noamb, aes(x=ChromosomeState, y=f_wCpG_OddsRatioGCtoAT))+geom_boxplot(outlier.shape=NA)+ylab("Odds Ratio: GC to AT")+xlab("Chromosome: Ancestral | Current") + geom_hline(yintercept=1, linetype="dotted", color="grey", size=1)
f_wCpG_ylim1=boxplot.stats(Isochores100KB_wCpG.noamb$f_wCpG_OddsRatioGCtoAT)$stats[c(1, 5)]
f_wCpG_BoxPlot=f_wCpG_BoxPlot+coord_flip(ylim = f_wCpG_ylim1*1.1)
f_wCpG_BoxPlot
ggsave("Fig6C_f_wCpG_BoxPlot.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Intermediate", scale = 1, width = 3.5, height = 3, units = c("in"), dpi =1200, limitsize = TRUE)

f_nCpG_BoxPlot=ggplot(Isochores100KB_nCpG.noamb, aes(x=ChromosomeState, y=f_nCpG_OddsRatioGCtoAT))+geom_boxplot(outlier.shape=NA)+ylab("Odds Ratio: GC to AT")+xlab("Chromosome: Ancestral | Current") + geom_hline(yintercept=1, linetype="dotted", color="grey", size=1)
f_nCpG_ylim1=boxplot.stats(Isochores100KB_nCpG.noamb$f_nCpG_OddsRatioGCtoAT)$stats[c(1, 5)]
f_nCpG_BoxPlot=f_nCpG_BoxPlot+coord_flip(ylim = f_nCpG_ylim1*1.1)
f_nCpG_BoxPlot
ggsave("Fig6D_f_nCpG_BoxPlot.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Intermediate", scale = 1, width = 3.5, height = 3, units = c("in"), dpi =1200, limitsize = TRUE)


Isochores100KB$Indel_Net=Isochores100KB$InsertionBP-Isochores100KB$DeletionBP
Isochores100KB$Indel_Sum=Isochores100KB$InsertionBP+Isochores100KB$DeletionBP
Isochores100KB.noamb=Isochores100KB[Isochores100KB$FalcoCONSENSUS_CLASS_MATCH!="Ambiguous",]
Isochores100KB.noamb=Isochores100KB.noamb[Isochores100KB.noamb$CONSENSUS_CLASS_MATCH!="Ambiguous",]

Isochores100KB.noamb$ANCESTRAL_STATE=rep("NA", length(Isochores100KB.noamb$CONSENSUS_CLASS_MATCH))
Isochores100KB.noamb$ANCESTRAL_STATE[Isochores100KB.noamb$CONSENSUS_CLASS_MATCH!="Microchromosome"]="Large"
Isochores100KB.noamb$ANCESTRAL_STATE[Isochores100KB.noamb$CONSENSUS_CLASS_MATCH=="Microchromosome"]="Microchromosome"

Isochores100KB.noamb$FALCO_STATE=rep("NA", length(Isochores100KB.noamb$FalcoCONSENSUS_CLASS_MATCH))
Isochores100KB.noamb$FALCO_STATE[Isochores100KB.noamb$FalcoCONSENSUS_CLASS_MATCH!="Microchromosome"]="Large"
Isochores100KB.noamb$FALCO_STATE[Isochores100KB.noamb$FalcoCONSENSUS_CLASS_MATCH=="Microchromosome"]="Microchromosome"

Isochores100KB.noamb$ChromosomeState=paste(Isochores100KB.noamb$ANCESTRAL_STATE, Isochores100KB.noamb$FALCO_STATE, sep=" | ")
Isochores100KB.noamb$ChromosomeState=gsub("Microchromosome","Micro",Isochores100KB.noamb$ChromosomeState)
Isochores100KB.noamb$ChromosomeState=factor(Isochores100KB.noamb$ChromosomeState, c("Large | Micro", "Large | Large", "Micro | Large", "Micro | Micro"))

pGC.anova=aov(pGC~Genome+ANCESTRAL_STATE:FALCO_STATE, data=Isochores100KB.noamb)
pGC_Tukey=TukeyHSD(x=pGC.anova, 'ANCESTRAL_STATE:FALCO_STATE', conf.level=0.95)
summary(pGC.anova)
pGC_Tukey
png("suppFig8_pGC_Tukey.png", width=14, height=7, units="in", res=1200)
par(mar=c(5,26,4,1)+.1)
plot(pGC_Tukey, las=1 , col="brown", cex=0.5)
dev.off()
pGC_BoxPlot=ggplot(Isochores100KB.noamb, aes(x=ChromosomeState, y=pGC*100))+geom_boxplot(outlier.shape=NA)+ylab("Percent GC")+xlab("Chromosome: Ancestral | Current")
pGC_ylim1 = boxplot.stats(Isochores100KB.noamb$pGC*100)$stats[c(1, 5)]
pGC_BoxPlot=pGC_BoxPlot+coord_flip(ylim = c(30,pGC_ylim1*1.2))
pGC_BoxPlot
ggsave("Fig6A_pGC_BoxPlot.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Intermediate", scale = 1, width = 3.5, height = 3, units = c("in"), dpi =1200, limitsize = TRUE)

CpG_nogccontrol.anova=aov(log2(CpG/tBASES+1)~Genome+ANCESTRAL_STATE:FALCO_STATE, data=Isochores100KB.noamb)
CpG_nogccontrol_Tukey=TukeyHSD(x=CpG_nogccontrol.anova, 'ANCESTRAL_STATE:FALCO_STATE', conf.level=0.95)
summary(CpG_nogccontrol.anova)
CpG_nogccontrol_Tukey
png("suppFig10CpG_nogccontrol_Tukey.png", width=14, height=7, units="in", res=1200)
par(mar=c(5,26,4,1)+.1)
plot(CpG_nogccontrol_Tukey, las=1 , col="brown", cex=0.5)
dev.off()
CpG_BoxPlot=ggplot(Isochores100KB.noamb, aes(x=ChromosomeState, y=CpG/tBASES*100))+geom_boxplot(outlier.shape=NA)+ylab("Percent CpG Sites")+xlab("Chromosome: Ancestral | Current")
CpG_ylim1 = boxplot.stats(Isochores100KB.noamb$CpG/Isochores100KB.noamb$tBASES*100)$stats[c(1, 5)]
CpG_BoxPlot=CpG_BoxPlot+coord_flip(ylim =CpG_ylim1*9)
CpG_BoxPlot
ggsave("Fig6B_GpC_BoxPlot.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Intermediate", scale = 1, width = 3.5, height = 3, units = c("in"), dpi =1200, limitsize = TRUE)

CpG.anova.gc=aov(log2(CpG/tBASES+1)~pGC, data=Isochores100KB.noamb)
CpG.anova.gcres=residuals(CpG.anova.gc)
ad.test(CpG.anova.gcres)
CpG.anova=step(aov(CpG.anova.gcres~Genome+ANCESTRAL_STATE:FALCO_STATE, data=Isochores100KB.noamb), direction="both")
CpG.anova.resid=residuals(CpG.anova)
ad.test(CpG.anova.resid)
CpG_Tukey=TukeyHSD(x=CpG.anova, 'ANCESTRAL_STATE:FALCO_STATE', conf.level=0.95)
summary(CpG.anova.gc)
summary(CpG.anova)
CpG_Tukey
png("suppFig9_CpG_Tukey.png", width=14, height=7, units="in", res=1200)
par(mar=c(5,26,4,1)+.1)
plot(CpG_Tukey, las=1 , col="brown", cex=0.5)
dev.off()
coef(summary.lm(CpG.anova.gc))

NetIndel.anova=step(aov(Indel_Net~Genome+Family+ANCESTRAL_STATE:FALCO_STATE, data=Isochores100KB.noamb), direction="both")
NetIndel.anova.res=residuals(NetIndel.anova)
ad.test(NetIndel.anova.res)
NetIndel_Tukey=TukeyHSD(x=NetIndel.anova, 'ANCESTRAL_STATE:FALCO_STATE', conf.level=0.95)
summary(NetIndel.anova)
NetIndel_Tukey
png("suppFig12_NetIndel_Tukey.png", width=14, height=7, units="in", res=1200)
par(mar=c(5,26,4,1)+.1)
plot(NetIndel_Tukey, las=1 , col="brown", cex=0.5)
dev.off()
NetIndel_BoxPlot=ggplot(Isochores100KB.noamb, aes(x=ChromosomeState, y=Indel_Net))+geom_boxplot(outlier.shape=NA)+ylab("Net Inserted and Deleted Bases")+xlab("Chromosome: Ancestral | Current")
NetIndel_ylim1 = boxplot.stats(Isochores100KB.noamb$Indel_Net)$stats[c(1, 5)]
NetIndel_BoxPlot=NetIndel_BoxPlot+coord_flip(ylim=c(-6, 12))
NetIndel_BoxPlot
ggsave("Fig6E_NetIndel_BoxPlot.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Intermediate", scale = 1, width = 3.5, height = 3, units = c("in"), dpi =1200, limitsize = TRUE)

#Without Controlling for GC
Isochores100KB.noamb$RepBP=Isochores100KB.noamb$BPrDNA+Isochores100KB.noamb$BPrLTR+Isochores100KB.noamb$BPrLINE+Isochores100KB.noamb$BPrSINE+Isochores100KB.noamb$BPrUnknown
NetRepBP.nogc.anova=step(aov(log2(RepBP+1)~Genome+ANCESTRAL_STATE:FALCO_STATE, data=Isochores100KB.noamb[is.finite(Isochores100KB.noamb$RepBP),]), direction="both")
NetRepBP.nogc.anova.resid=residuals(NetRepBP.nogc.anova)
ad.test(NetRepBP.nogc.anova.resid)
NetRepBPnoGC_Tukey=TukeyHSD(x=NetRepBP.nogc.anova, 'ANCESTRAL_STATE:FALCO_STATE', conf.level=0.95)
summary(NetRepBP.nogc.anova)
NetRepBPnoGC_Tukey
png("suppFig13_NetRepBPnoGC_Tukey.png", width=14, height=7, units="in", res=1200)
par(mar=c(5,26,4,1)+.1)
plot(NetRepBPnoGC_Tukey, las=1 , col="brown", cex=0.5)
dev.off()
coef(summary.lm(NetRepBP.nogc.anova))
NetRepBP_BoxPlot=ggplot(Isochores100KB.noamb, aes(x=ChromosomeState, y=RepBP))+geom_boxplot(outlier.shape=NA)+ylab("Repetitive Element Bases")+xlab("Chromosome: Ancestral | Current")
NetRepBP_ylim1 = boxplot.stats(Isochores100KB.noamb$RepBP)$stats[c(1, 5)]
NetRepBP_BoxPlot=NetRepBP_BoxPlot+coord_flip(ylim=NetRepBP_ylim1*1.1)
NetRepBP_BoxPlot
ggsave("Fig6F_NetRepBP_BoxPlot.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Intermediate", scale = 1, width = 3.5, height = 3, units = c("in"), dpi =1200, limitsize = TRUE)

NetRepBP.anova.gc=aov(log2(RepBP+1)~pGC, data=Isochores100KB.noamb)
NetRepBP.anova.gc.resid=residuals(NetRepBP.anova.gc)
ad.test(NetRepBP.anova.gc.resid)
NetRepBP.anova=step(aov(NetRepBP.anova.gc.resid~Genome+ANCESTRAL_STATE:FALCO_STATE, data=Isochores100KB.noamb), direction="both")
NetRepBP.anova.resid=residuals(NetRepBP.anova)
ad.test(NetRepBP.anova.resid)
NetRepBP_Tukey=TukeyHSD(x=NetRepBP.anova, 'ANCESTRAL_STATE:FALCO_STATE', conf.level=0.95)
summary(NetRepBP.anova)
NetRepBP_Tukey
png("suppFig14_NetRepBP_Tukey.png", width=14, height=7, units="in", res=1200)
par(mar=c(5,26,4,1)+.1)
plot(NetRepBP_Tukey, las=1 , col="brown", cex=0.5)
dev.off()
coef(summary.lm(NetRepBP.anova.gc))


setwd("~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Supplementary/Unused")

SumIndel.anova.gc=step(aov(log2(Indel_Sum+1)~pGC, data=Isochores100KB.noamb), direction="both")
SumIndel.anova.gcres=residuals(SumIndel.anova.gc)
ad.test(SumIndel.anova.gcres)
Isochores100KB.noamb.indel=Isochores100KB.noamb[row.names(as.data.frame(SumIndel.anova.gcres)),]
SumIndel.anova=step(aov(SumIndel.anova.gcres~Genome+ANCESTRAL_STATE:FALCO_STATE, data=Isochores100KB.noamb.indel), direction="both")
SumIndel.anova.res=residuals(SumIndel.anova)
SumIndel_Tukey=TukeyHSD(x=SumIndel.anova, 'ANCESTRAL_STATE:FALCO_STATE', conf.level=0.95)
summary(SumIndel.anova)
SumIndel_Tukey
png("Unused_SumIndel_Tukey.png", width=14, height=7, units="in", res=1200)
par(mar=c(5,26,4,1)+.1)
plot(SumIndel_Tukey, las=1 , col="brown", cex=0.5)
dev.off()
coef(summary.lm(SumIndel.anova.gc))

#Sex Chromosome
Isochores100KB.noamb_sex=Isochores100KB[Isochores100KB$ConFalcoSex_Chrom_Type!="Ambiguous",]
Isochores100KB.noamb_sex$RepBP=Isochores100KB.noamb_sex$BPrDNA+Isochores100KB.noamb_sex$BPrLTR+Isochores100KB.noamb_sex$BPrLINE+Isochores100KB.noamb_sex$BPrSINE+Isochores100KB.noamb_sex$BPrUnknown
SexRepeats.anova.gc=aov(log2(RepBP+1)~pGC, data=Isochores100KB.noamb_sex)
SexRepeats.anova.gc.resid=residuals(SexRepeats.anova.gc)
ad.test(SexRepeats.anova.gc.resid)
SexRepeats.anova=aov(SexRepeats.anova.gc.resid~ConFalcoSex_Chrom_Type, data=Isochores100KB.noamb_sex)
SexRepeats.anova.resid=residuals(SexRepeats.anova)
ad.test(SexRepeats.anova.resid)
SexRepeats_Tukey=TukeyHSD(x=SexRepeats.anova, 'ConFalcoSex_Chrom_Type', conf.level=0.95)
SexRepeats_Tukey
png("Unused_SexRepeats_Tukey.png", width=14, height=7, units="in", res=1200)
par(mar=c(5,26,4,1)+.1)
plot(SexRepeats_Tukey, las=1 , col="brown", cex=0.5)
dev.off()
coef(summary.lm(SexRepeats.anova))
summary(SexRepeats.anova)


Isochores100KB.noamb$SV_Sum=Isochores100KB.noamb$BPuInsSV+Isochores100KB.noamb$BPuDelSV+Isochores100KB.noamb$BPuInvSV+Isochores100KB.noamb$BPuSeg_DelSV+Isochores100KB.noamb$BPuSeg_InsSV
SV_Sum.anova.gc=aov(log2(SV_Sum+1)~pGC, data=Isochores100KB.noamb)
SV_Sum.anova.gcres=residuals(SV_Sum.anova.gc)
ad.test(SV_Sum.anova.gcres)
Isochores100KB.noamb.sv=Isochores100KB.noamb[row.names(as.data.frame(SV_Sum.anova.gcres)),]
SV_Sum.anova=aov(SV_Sum.anova.gcres~Genome+ANCESTRAL_STATE:FALCO_STATE, data=Isochores100KB.noamb.sv)
SV_Sum.anova.res=residuals(SV_Sum.anova)
ad.test(SV_Sum.anova.res)
SV_Sum_Tukey=TukeyHSD(x=SV_Sum.anova, 'ANCESTRAL_STATE:FALCO_STATE', conf.level=0.95)
summary(SV_Sum.anova)
SV_Sum_Tukey
png("Unused_SV_Sum_Tukey.png", width=14, height=7, units="in", res=1200)
par(mar=c(5,26,4,1)+.1)
plot(SV_Sum_Tukey, las=1 , col="brown", cex=0.5)
dev.off()
coef(summary.lm(SV_Sum.anova.gc))

Isochores100KB.noamb$SV_Net=Isochores100KB.noamb$BPuInsSV-Isochores100KB.noamb$BPuDelSV+Isochores100KB.noamb$BPuInvSV-Isochores100KB.noamb$BPuSeg_DelSV+Isochores100KB.noamb$BPuSeg_InsSV
SV_Net.anova.gc=aov(SV_Net~pGC, data=Isochores100KB.noamb)
SV_Net.anova.gcres=residuals(SV_Net.anova.gc)
ad.test(SV_Net.anova.gcres)
Isochores100KB.noamb.sv=Isochores100KB.noamb[row.names(as.data.frame(SV_Net.anova.gcres)),]
SV_Net.anova=aov(SV_Net.anova.gcres~Genome+ANCESTRAL_STATE:FALCO_STATE, data=Isochores100KB.noamb.sv)
SV_Net.anova.res=residuals(SV_Net.anova)
ad.test(SV_Net.anova.res)
SV_Net_Tukey=TukeyHSD(x=SV_Net.anova, 'ANCESTRAL_STATE:FALCO_STATE', conf.level=0.95)
summary(SV_Net.anova)
SV_Net_Tukey
png("Unused_SV_Net_Tukey.png", width=14, height=7, units="in", res=1200)
par(mar=c(5,26,4,1)+.1)
plot(SV_Net_Tukey, las=1 , col="brown", cex=0.5)
dev.off()
coef(summary.lm(SV_Net.anova.gc))

ROH.anova.gc=aov(ROH~pGC, data=Isochores100KB.noamb)
ROH.anova=step(aov(ROH~Genome+ANCESTRAL_STATE:FALCO_STATE+Family, data=Isochores100KB.noamb), direction="both")
ROH_Tukey=TukeyHSD(x=ROH.anova, 'ANCESTRAL_STATE:FALCO_STATE', conf.level=0.95)
summary(ROH.anova)
ROH_Tukey
png("Unused_Tukey.png", width=14, height=7, units="in", res=1200)
par(mar=c(5,26,4,1)+.1)
plot(ROH_Tukey, las=1 , col="brown", cex=0.5)
dev.off()
coef(summary.lm(ROH.anova.gc))
#######Create Table of MUMmer Alignments to Chromosome Types#########################
setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis/Synthesis/Data/ChromAnnotated")
Isochores100KB=read.delim("AllultChromAnnotated100KBWindowsSummary.txt", header=TRUE, sep="\t")
Isochores100KB$Genome=as.character(Isochores100KB$Genome)
Isochores100KB$Genome[Isochores100KB$Genome=="A"]="Lanner"
Isochores100KB$Genome[Isochores100KB$Genome=="B"]="Barbary"
Isochores100KB$Genome[Isochores100KB$Genome=="C"]="Peregrine"
Isochores100KB$Genome[Isochores100KB$Genome=="D"]="Saker"
Isochores100KB$Genome[Isochores100KB$Genome=="E"]="Gyr-1"
Isochores100KB$Genome[Isochores100KB$Genome=="F"]="Gyr-2"
Isochores100KB$Genome[Isochores100KB$Genome=="FpBSx77"]="Bl. Shaheen"
Isochores100KB$Genome[Isochores100KB$Genome=="Kestrel"]="C. Kestrel"


Lanner_Isochores100KB=Isochores100KB[Isochores100KB$Genome=="Lanner",]
Barbary_Isochores100KB=Isochores100KB[Isochores100KB$Genome=="Barbary",]
Peregrine_Isochores100KB=Isochores100KB[Isochores100KB$Genome=="Peregrine",]
Saker_Isochores100KB=Isochores100KB[Isochores100KB$Genome=="Saker",]
Gyr1_Isochores100KB=Isochores100KB[Isochores100KB$Genome=="Gyr-1",]
Gyr2_Isochores100KB=Isochores100KB[Isochores100KB$Genome=="Gyr-2",]
BlackShaheen_Isochores100KB=Isochores100KB[Isochores100KB$Genome=="Bl. Shaheen",]
Kestrel_Isochores100KB=Isochores100KB[Isochores100KB$Genome=="C. Kestrel",]

LannerChrom=data.frame(Genome="Lanner")
BarbaryChrom=data.frame(Genome="Barbary")
PeregrineChrom=data.frame(Genome="Peregrine")
SakerChrom=data.frame(Genome="Saker")
Gyr1Chrom=data.frame(Genome="Gyr-1")
Gyr2Chrom=data.frame(Genome="Gyr-2")
BlackShaheenChrom=data.frame(Genome="Bl. Shaheen")
KestrelChrom=data.frame(Genome="C. Kestrel")

LannerChrom$Windows=dim(Lanner_Isochores100KB)[1]
BarbaryChrom$Windows=dim(Barbary_Isochores100KB)[1]
PeregrineChrom$Windows=dim(Peregrine_Isochores100KB)[1]
SakerChrom$Windows=dim(Saker_Isochores100KB)[1]
Gyr1Chrom$Windows=dim(Gyr1_Isochores100KB)[1]
Gyr2Chrom$Windows=dim(Gyr2_Isochores100KB)[1]
BlackShaheenChrom$Windows=dim(BlackShaheen_Isochores100KB)[1]
KestrelChrom$Windows=dim(Kestrel_Isochores100KB)[1]

LannerChrom$FalcoMacro=dim(Lanner_Isochores100KB[Lanner_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
BarbaryChrom$FalcoMacro=dim(Barbary_Isochores100KB[Barbary_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
PeregrineChrom$FalcoMacro=dim(Peregrine_Isochores100KB[Peregrine_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
SakerChrom$FalcoMacro=dim(Saker_Isochores100KB[Saker_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
Gyr1Chrom$FalcoMacro=dim(Gyr1_Isochores100KB[Gyr1_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
Gyr2Chrom$FalcoMacro=dim(Gyr2_Isochores100KB[Gyr2_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
BlackShaheenChrom$FalcoMacro=dim(BlackShaheen_Isochores100KB[BlackShaheen_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
KestrelChrom$FalcoMacro=dim(Kestrel_Isochores100KB[Kestrel_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]

LannerChrom$FalcoMacro=dim(Lanner_Isochores100KB[Lanner_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
BarbaryChrom$FalcoMacro=dim(Barbary_Isochores100KB[Barbary_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
PeregrineChrom$FalcoMacro=dim(Peregrine_Isochores100KB[Peregrine_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
SakerChrom$FalcoMacro=dim(Saker_Isochores100KB[Saker_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
Gyr1Chrom$FalcoMacro=dim(Gyr1_Isochores100KB[Gyr1_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
Gyr2Chrom$FalcoMacro=dim(Gyr2_Isochores100KB[Gyr2_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
BlackShaheenChrom$FalcoMacro=dim(BlackShaheen_Isochores100KB[BlackShaheen_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
KestrelChrom$FalcoMacro=dim(Kestrel_Isochores100KB[Kestrel_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]

LannerChrom$FalcoIntermediate=dim(Lanner_Isochores100KB[Lanner_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Intermediate",])[1]
BarbaryChrom$FalcoIntermediate=dim(Barbary_Isochores100KB[Barbary_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Intermediate",])[1]
PeregrineChrom$FalcoIntermediate=dim(Peregrine_Isochores100KB[Peregrine_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Intermediate",])[1]
SakerChrom$FalcoIntermediate=dim(Saker_Isochores100KB[Saker_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Intermediate",])[1]
Gyr1Chrom$FalcoIntermediate=dim(Gyr1_Isochores100KB[Gyr1_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Intermediate",])[1]
Gyr2Chrom$FalcoIntermediate=dim(Gyr2_Isochores100KB[Gyr2_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Intermediate",])[1]
BlackShaheenChrom$FalcoIntermediate=dim(BlackShaheen_Isochores100KB[BlackShaheen_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Intermediate",])[1]
KestrelChrom$FalcoIntermediate=dim(Kestrel_Isochores100KB[Kestrel_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Intermediate",])[1]

LannerChrom$FalcoMicro=dim(Lanner_Isochores100KB[Lanner_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
BarbaryChrom$FalcoMicro=dim(Barbary_Isochores100KB[Barbary_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
PeregrineChrom$FalcoMicro=dim(Peregrine_Isochores100KB[Peregrine_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
SakerChrom$FalcoMicro=dim(Saker_Isochores100KB[Saker_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
Gyr1Chrom$FalcoMicro=dim(Gyr1_Isochores100KB[Gyr1_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
Gyr2Chrom$FalcoMicro=dim(Gyr2_Isochores100KB[Gyr2_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
BlackShaheenChrom$FalcoMicro=dim(BlackShaheen_Isochores100KB[BlackShaheen_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
KestrelChrom$FalcoMicro=dim(Kestrel_Isochores100KB[Kestrel_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Microchromosome",])[1]

LannerChrom$FalcoAmbiguous=dim(Lanner_Isochores100KB[Lanner_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Ambiguous",])[1]
BarbaryChrom$FalcoAmbiguous=dim(Barbary_Isochores100KB[Barbary_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Ambiguous",])[1]
PeregrineChrom$FalcoAmbiguous=dim(Peregrine_Isochores100KB[Peregrine_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Ambiguous",])[1]
SakerChrom$FalcoAmbiguous=dim(Saker_Isochores100KB[Saker_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Ambiguous",])[1]
Gyr1Chrom$FalcoAmbiguous=dim(Gyr1_Isochores100KB[Gyr1_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Ambiguous",])[1]
Gyr2Chrom$FalcoAmbiguous=dim(Gyr2_Isochores100KB[Gyr2_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Ambiguous",])[1]
BlackShaheenChrom$FalcoAmbiguous=dim(BlackShaheen_Isochores100KB[BlackShaheen_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Ambiguous",])[1]
KestrelChrom$FalcoAmbiguous=dim(Kestrel_Isochores100KB[Kestrel_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="Ambiguous",])[1]

LannerChrom$FalcoUnaligned=dim(Lanner_Isochores100KB[Lanner_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="",])[1]
BarbaryChrom$FalcoUnaligned=dim(Barbary_Isochores100KB[Barbary_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="",])[1]
PeregrineChrom$FalcoUnaligned=dim(Peregrine_Isochores100KB[Peregrine_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="",])[1]
SakerChrom$FalcoUnaligned=dim(Saker_Isochores100KB[Saker_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="",])[1]
Gyr1Chrom$FalcoUnaligned=dim(Gyr1_Isochores100KB[Gyr1_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="",])[1]
Gyr2Chrom$FalcoUnaligned=dim(Gyr2_Isochores100KB[Gyr2_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="",])[1]
BlackShaheenChrom$FalcoUnaligned=dim(BlackShaheen_Isochores100KB[BlackShaheen_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="",])[1]
KestrelChrom$FalcoUnaligned=dim(Kestrel_Isochores100KB[Kestrel_Isochores100KB$FalcoCONSENSUS_CLASS_MATCH=="",])[1]


LannerChrom$AncestralMacro=dim(Lanner_Isochores100KB[Lanner_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
BarbaryChrom$AncestralMacro=dim(Barbary_Isochores100KB[Barbary_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
PeregrineChrom$AncestralMacro=dim(Peregrine_Isochores100KB[Peregrine_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
SakerChrom$AncestralMacro=dim(Saker_Isochores100KB[Saker_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
Gyr1Chrom$AncestralMacro=dim(Gyr1_Isochores100KB[Gyr1_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
Gyr2Chrom$AncestralMacro=dim(Gyr2_Isochores100KB[Gyr2_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
BlackShaheenChrom$AncestralMacro=dim(BlackShaheen_Isochores100KB[BlackShaheen_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
KestrelChrom$AncestralMacro=dim(Kestrel_Isochores100KB[Kestrel_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]

LannerChrom$AncestralMacro=dim(Lanner_Isochores100KB[Lanner_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
BarbaryChrom$AncestralMacro=dim(Barbary_Isochores100KB[Barbary_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
PeregrineChrom$AncestralMacro=dim(Peregrine_Isochores100KB[Peregrine_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
SakerChrom$AncestralMacro=dim(Saker_Isochores100KB[Saker_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
Gyr1Chrom$AncestralMacro=dim(Gyr1_Isochores100KB[Gyr1_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
Gyr2Chrom$AncestralMacro=dim(Gyr2_Isochores100KB[Gyr2_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
BlackShaheenChrom$AncestralMacro=dim(BlackShaheen_Isochores100KB[BlackShaheen_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]
KestrelChrom$AncestralMacro=dim(Kestrel_Isochores100KB[Kestrel_Isochores100KB$CONSENSUS_CLASS_MATCH=="Macrochromosome",])[1]

LannerChrom$AncestralIntermediate=dim(Lanner_Isochores100KB[Lanner_Isochores100KB$CONSENSUS_CLASS_MATCH=="Intermediate",])[1]
BarbaryChrom$AncestralIntermediate=dim(Barbary_Isochores100KB[Barbary_Isochores100KB$CONSENSUS_CLASS_MATCH=="Intermediate",])[1]
PeregrineChrom$AncestralIntermediate=dim(Peregrine_Isochores100KB[Peregrine_Isochores100KB$CONSENSUS_CLASS_MATCH=="Intermediate",])[1]
SakerChrom$AncestralIntermediate=dim(Saker_Isochores100KB[Saker_Isochores100KB$CONSENSUS_CLASS_MATCH=="Intermediate",])[1]
Gyr1Chrom$AncestralIntermediate=dim(Gyr1_Isochores100KB[Gyr1_Isochores100KB$CONSENSUS_CLASS_MATCH=="Intermediate",])[1]
Gyr2Chrom$AncestralIntermediate=dim(Gyr2_Isochores100KB[Gyr2_Isochores100KB$CONSENSUS_CLASS_MATCH=="Intermediate",])[1]
BlackShaheenChrom$AncestralIntermediate=dim(BlackShaheen_Isochores100KB[BlackShaheen_Isochores100KB$CONSENSUS_CLASS_MATCH=="Intermediate",])[1]
KestrelChrom$AncestralIntermediate=dim(Kestrel_Isochores100KB[Kestrel_Isochores100KB$CONSENSUS_CLASS_MATCH=="Intermediate",])[1]

LannerChrom$AncestralMicro=dim(Lanner_Isochores100KB[Lanner_Isochores100KB$CONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
BarbaryChrom$AncestralMicro=dim(Barbary_Isochores100KB[Barbary_Isochores100KB$CONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
PeregrineChrom$AncestralMicro=dim(Peregrine_Isochores100KB[Peregrine_Isochores100KB$CONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
SakerChrom$AncestralMicro=dim(Saker_Isochores100KB[Saker_Isochores100KB$CONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
Gyr1Chrom$AncestralMicro=dim(Gyr1_Isochores100KB[Gyr1_Isochores100KB$CONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
Gyr2Chrom$AncestralMicro=dim(Gyr2_Isochores100KB[Gyr2_Isochores100KB$CONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
BlackShaheenChrom$AncestralMicro=dim(BlackShaheen_Isochores100KB[BlackShaheen_Isochores100KB$CONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
KestrelChrom$AncestralMicro=dim(Kestrel_Isochores100KB[Kestrel_Isochores100KB$CONSENSUS_CLASS_MATCH=="Microchromosome",])[1]

LannerChrom$AncestralAmbiguous=dim(Lanner_Isochores100KB[Lanner_Isochores100KB$CONSENSUS_CLASS_MATCH=="Ambiguous",])[1]
BarbaryChrom$AncestralAmbiguous=dim(Barbary_Isochores100KB[Barbary_Isochores100KB$CONSENSUS_CLASS_MATCH=="Ambiguous",])[1]
PeregrineChrom$AncestralAmbiguous=dim(Peregrine_Isochores100KB[Peregrine_Isochores100KB$CONSENSUS_CLASS_MATCH=="Ambiguous",])[1]
SakerChrom$AncestralAmbiguous=dim(Saker_Isochores100KB[Saker_Isochores100KB$CONSENSUS_CLASS_MATCH=="Ambiguous",])[1]
Gyr1Chrom$AncestralAmbiguous=dim(Gyr1_Isochores100KB[Gyr1_Isochores100KB$CONSENSUS_CLASS_MATCH=="Ambiguous",])[1]
Gyr2Chrom$AncestralAmbiguous=dim(Gyr2_Isochores100KB[Gyr2_Isochores100KB$CONSENSUS_CLASS_MATCH=="Ambiguous",])[1]
BlackShaheenChrom$AncestralAmbiguous=dim(BlackShaheen_Isochores100KB[BlackShaheen_Isochores100KB$CONSENSUS_CLASS_MATCH=="Ambiguous",])[1]
KestrelChrom$AncestralAmbiguous=dim(Kestrel_Isochores100KB[Kestrel_Isochores100KB$CONSENSUS_CLASS_MATCH=="Ambiguous",])[1]

LannerChrom$AncestralUnaligned=dim(Lanner_Isochores100KB[Lanner_Isochores100KB$CONSENSUS_CLASS_MATCH=="",])[1]
BarbaryChrom$AncestralUnaligned=dim(Barbary_Isochores100KB[Barbary_Isochores100KB$CONSENSUS_CLASS_MATCH=="",])[1]
PeregrineChrom$AncestralUnaligned=dim(Peregrine_Isochores100KB[Peregrine_Isochores100KB$CONSENSUS_CLASS_MATCH=="",])[1]
SakerChrom$AncestralUnaligned=dim(Saker_Isochores100KB[Saker_Isochores100KB$CONSENSUS_CLASS_MATCH=="",])[1]
Gyr1Chrom$AncestralUnaligned=dim(Gyr1_Isochores100KB[Gyr1_Isochores100KB$CONSENSUS_CLASS_MATCH=="",])[1]
Gyr2Chrom$AncestralUnaligned=dim(Gyr2_Isochores100KB[Gyr2_Isochores100KB$CONSENSUS_CLASS_MATCH=="",])[1]
BlackShaheenChrom$AncestralUnaligned=dim(BlackShaheen_Isochores100KB[BlackShaheen_Isochores100KB$CONSENSUS_CLASS_MATCH=="",])[1]
KestrelChrom$AncestralUnaligned=dim(Kestrel_Isochores100KB[Kestrel_Isochores100KB$CONSENSUS_CLASS_MATCH=="",])[1]

AllFalco100KBWindowChromAlign=rbind(LannerChrom,SakerChrom,PeregrineChrom,BarbaryChrom,BlackShaheenChrom,Gyr1Chrom,Gyr2Chrom,KestrelChrom)

AllFalcoPercentAlign=data.frame(Genome=AllFalco100KBWindowChromAlign$Genome, Windows=AllFalco100KBWindowChromAlign$Windows, FalcoUnaligned=AllFalco100KBWindowChromAlign$FalcoUnaligned, AncestralUnaligned=AllFalco100KBWindowChromAlign$AncestralUnaligned)
AllFalcoPercentAlign$pFalcoUnaligned=AllFalcoPercentAlign$FalcoUnaligned/AllFalcoPercentAlign$Windows*100
AllFalcoPercentAlign$pAncestralUnaligned=AllFalcoPercentAlign$AncestralUnaligned/AllFalcoPercentAlign$Windows*100
AllFalcoPercentAlign$pAncestralAligned=100-AllFalcoPercentAlign$pAncestralUnaligned

setwd("~/NYU_AD/FalconProject/Manuscripts/FalconRef8/UltimateTables/SupplementaryTables/Intermediate")
write.table(AllFalco100KBWindowChromAlign, file="SupplementaryTable_100KBChromTypeAlign.txt", sep="\t", row.names=FALSE)
setwd("~/NYU_AD/FalconProject/Manuscripts/FalconRef8/UltimateTables/Intermediate")
write.table(AllFalcoPercentAlign, file="ChromAlignmentPercents.txt", sep="\t", row.names=FALSE)
#############Create NUMT Figues##################
#Running Analysis of NUMTs on pennultimate clusters
#Statistical Analysis of NUMTs in R Studio
#Import trees and tables
setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis/ultMitochondrial/NUMTs/Trees/ultClusters/UltimateTrees")
library(ape)
for(i in 0:46){
assign(gsub(" ", "", paste("Cluster",i,".tre")), read.tree(gsub(" ", "", paste("RAxML_bipartitionsBranchLabels.Cluster",i,".tre"))))
}

setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis/ultMitochondrial/NUMTs/Trees/ultClusters/UltimateTables")
for(i in 0:46){
assign(gsub(" ", "", paste("Cluster",i,"_Table")), read.delim(gsub(" ", "", paste("RAxML_bipartitionsBranchLabels.Cluster",i,"_Table",".txt")), sep="\t", row.names=1))
}

#Extract Tree Distances with APE

#Create vectors of Tree and Table Objects
TreeNames=c(rep(NA, 47))
TableNames=c(rep(NA, 47))
Count=1
for (i in 0:46){
TreeNames[Count]=gsub(" ", "", paste("Cluster",i,".tre"))
TableNames[Count]=gsub(" ", "", paste("Cluster",i,"_Table"))
Count=Count+1
}

for (i in 0:46){
v=TreeNames[i+1]
assign(gsub(" ", "", paste("Cluster",i,".","dist")), as.data.frame(cophenetic.phylo(get(v))))
}

#Extract Table Information: Basepair Length Distances
library(vegan)
for (i in 0:46){
t=TableNames[i+1]
assign(gsub(" ", "", paste("Cluster",i,"basepairs.euclid")), vegdist(get(t)$BASEPAIRS[!is.na(get(t)$BASEPAIRS)], method="euclid"))
}

for (i in 0:46){
t=TableNames[i+1]
assign(gsub(" ", "", paste("Cluster",i,"basepairs.euclid")), as.matrix(vegdist(get(t)$BASEPAIRS[!is.na(get(t)$BASEPAIRS)], method="euclidean")))
}


#Extract Variation in Trees Explained by Synteny
SyntenyTreeInfo=data.frame(matrix(nrow=0, ncol=5))
colnames(SyntenyTreeInfo)=c("Cluster","R2","NumberNumts", "p-value", "F-value")
for (i in 0:46){
v=TreeNames[i+1]
t=TableNames[i+1]
tree=get(gsub(" ", "", paste("Cluster",i,".","dist")))
synt=as.factor(get(t)$SYNT)
tree.ord=tree[row.names(get(t)),row.names(get(t))]
synt.purged=synt[!is.na(synt)]
tree.synt.ord.purged=tree.ord[synt[!is.na(synt)], synt[!is.na(synt)]]
TreeSynteny.adonis2=adonis2(tree.synt.ord.purged~synt.purged, nperm=1000)
TempSyntenyTreeInfo=cbind(i,TreeSynteny.adonis2$R2[1],length(synt.purged), TreeSynteny.adonis2$`Pr(>F)`[1], TreeSynteny.adonis2$F[1])
colnames(TempSyntenyTreeInfo)=c("Cluster","R2","NumberNumts", "p-value", "F-value")
SyntenyTreeInfo=rbind(SyntenyTreeInfo, TempSyntenyTreeInfo)
}
mean(SyntenyTreeInfo$R2[!is.na(SyntenyTreeInfo$R2)])
sd(SyntenyTreeInfo$R2[!is.na(SyntenyTreeInfo$R2)])
SyntenyTreeInfo$Weights=SyntenyTreeInfo$NumberNumts/sum(SyntenyTreeInfo$NumberNumts)
wmTreeSYNT=weighted.mean(SyntenyTreeInfo$R2[!is.na(SyntenyTreeInfo$R2)],SyntenyTreeInfo$Weights[!is.na(SyntenyTreeInfo$R2)])
sdTreeSYNT=sum(SyntenyTreeInfo$Weights[!is.na(SyntenyTreeInfo$R2)]*(SyntenyTreeInfo$R2[!is.na(SyntenyTreeInfo$R2)]-wmTreeSYNT)^2)
wmTreeSYNT
sdTreeSYNT

#Length Synteny Appendix
BASEPAIRSyntenyInfo=data.frame(matrix(nrow=0, ncol=5))
colnames(BASEPAIRSyntenyInfo)=c("Cluster","R2","NumberNumts", "p-value", "F-value")
for (i in 0:46){
t=TableNames[i+1]
AdonisTable=get(t)[!is.na(get(t)$BASEPAIRS),]
AdonisTable.purged=AdonisTable[!is.na(AdonisTable$SYNT),]
BASEPAIRSynteny.adonis2=adonis2(AdonisTable.purged$BASEPAIRS~AdonisTable.purged$SYNT,method="euclidean", nperm=1000)
TempBASEPAIRSyntenyInfo=cbind(i,BASEPAIRSynteny.adonis2$R2[1],length(AdonisTable.purged$SYNT), BASEPAIRSynteny.adonis2$`Pr(>F)`[1], BASEPAIRSynteny.adonis2$F[1])
colnames(TempBASEPAIRSyntenyInfo)=c("Cluster","R2","NumberNumts", "p-value", "F-value")
BASEPAIRSyntenyInfo=rbind(BASEPAIRSyntenyInfo, TempBASEPAIRSyntenyInfo)
}
mean(BASEPAIRSyntenyInfo$R2[!is.na(BASEPAIRSyntenyInfo$R2)])
sd(BASEPAIRSyntenyInfo$R2[!is.na(BASEPAIRSyntenyInfo$R2)])
BASEPAIRSyntenyInfo$Weights=BASEPAIRSyntenyInfo$NumberNumts/sum(BASEPAIRSyntenyInfo$NumberNumts)
wmBPSYNT=weighted.mean(BASEPAIRSyntenyInfo$R2[!is.na(BASEPAIRSyntenyInfo$R2)],BASEPAIRSyntenyInfo$Weights[!is.na(BASEPAIRSyntenyInfo$R2)])
sdBPSYNT=sum(BASEPAIRSyntenyInfo$Weights[!is.na(BASEPAIRSyntenyInfo$R2)]*(BASEPAIRSyntenyInfo$R2[!is.na(BASEPAIRSyntenyInfo$R2)]-wmBPSYNT)^2)
wmBPSYNT
sdBPSYNT

#Nearest Mitogenome
MergedNUMTtable=as.data.frame(matrix(nrow=0, ncol=14))
for (i in 0:46){
v=TreeNames[i+1]
t=TableNames[i+1]
tree=as.data.frame(get(gsub(" ", "", paste("Cluster",i,".","dist"))))
table=get(t)
mito=tree[,which(names(tree) %in% row.names(table[table$TYPE=="MITOGENOME",]))]
table$Cluster=i
table$nMITO=NA
table$sMitoDist=NA
for (n in row.names(tree)){
nearestMITO=names(which.min(mito[n,]))
self=gsub("([a-zA-Z]+)[0-9_]+.*","\\1Mitogenome",n)
self=gsub("(Gyr-1)[0-9_]+.*","\\1Mitogenome",self)
self=gsub("(Gyr-2)[0-9_]+.*","\\1Mitogenome",self)
selfMitoDist=mito[n,self]
table[n,]$nMITO=nearestMITO
table[n,]$sMitoDist=selfMitoDist
}
MergedNUMTtable=rbind(MergedNUMTtable, table)
}

MergedNUMTtable$Clade=NA
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="KestrelMitogenome"]="Kestrel"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="MW266990.1_Falco_subbuteo_mitochondrion"]="Mid-Falcon"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="MN431594.1_Falco_subbuteo_mitochondrion"]="Mid-Falcon"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_057089.1_Falco_subbuteo_mitochondrion"]="Mid-Falcon"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_026715.1_Falco_cherrug_mitochondrion"]="Hierofalco"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_029359.1_Falco_rusticolus_mitochondrion"]="Hierofalco"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="KT989235.1_Falco_rusticolus_mitochondrion"]="Hierofalco"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="KP337902.1_Falco_cherrug_mitochondrion"]="Hierofalco"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="KP337902.1_Falco_cherrug_mitochondrion"]="Hierofalco"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="LannerMitogenome"]="Hierofalco"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="Gyr-1Mitogenome"]="Hierofalco"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="Gyr-2Mitogenome"]="Hierofalco"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="SakerMitogenome"]="Hierofalco"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="CM026719.1_Falco_rusticolus_isolate_bFalRus1_mitochondrion"]="Hierofalco"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="JQ282801.1_Falco_peregrinus_mitochondrion"]="Peregrine"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="AF090338.1_Falco_peregrinus_mitochondrion"]="Peregrine"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="AF090338.1_Falco_peregrinus_mitochondrion"]="Peregrine"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_000878.1_Falco_peregrinus_mitochondrion"]="Peregrine"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="JX029991.1_Falco_peregrinus_mitochondrion"]="Peregrine"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="PeregrineMitogenome"]="Peregrine"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="BarbaryMitogenome"]="Peregrine"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="Bl.Shaheen"]="Peregrine"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="KX987839.1_Falco_amurensis_mitochondrion"]="Mid-Falcon"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_039842.1_Falco_amurensis_mitochondrion"]="Mid-Falcon"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="CM030192.1_Falco_naumanni_isolate_bFalNau1_mitochondrion"]="Kestrel"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_054728.1_Falco_naumanni_isolate_bFalNau1_mitochondrion"]="Kestrel"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="KM264304.1_Falco_columbarius_mitochondrion"]="Mid-Falcon"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_025579.1_Falco_columbarius_mitochondrion"]="Mid-Falcon"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="EU196361.1_Falco_tinnunculus_mitochondrion"]="Kestrel"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_011307.1_Falco_tinnunculus_mitochondrion"]="Kestrel"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="KM251414.1_Falco_naumanni_mitochondrion"]="Kestrel"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="DQ780880.1_Falco_sparverius_mitochondrion"]="Mid-Falcon"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_008547.1_Falco_sparverius_mitochondrion"]="Mid-Falcon"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_044674.1_Caracara_creightoni_voucher"]="Polyborinae"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="KP064202.1_Phalcoboenus_australis_mitochondrion"]="Polyborinae"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_031897.1_Phalcoboenus_australis_mitochondrion"]="Polyborinae"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="MN231450.1_Caracara_plancus_voucher"]="Polyborinae"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_044672.1_Caracara_plancus_voucher"]="Polyborinae"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_044672.1_Caracara_plancus_voucher"]="Polyborinae"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="MN231451.1_Caracara_cheriway_voucher"]="Polyborinae"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_044673.1_Caracara_cheriway_voucher"]="Polyborinae"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="MN231452.1_Caracara_creightoni_voucher"]="Polyborinae"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="DQ780881.1_Micrastur_gilvicollis_mitochondrion"]="Herpetherinae"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_008548.1_Micrastur_gilvicollis_mitochondrion"]="Herpetherinae"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="MN356294.1_Herpetotheres_cachinnans_mitochondrion"]="Herpetherinae"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_052801.1_Herpetotheres_cachinnans_mitochondrion"]="Herpetherinae"

MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="AP010797.1_Accipiter_gentilis_mitochondrial_DNA"]="Neoaves"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="CM020378.1_Catharus_ustulatus_isolate_bCatUst1_mitochondrion"]="Neoaves"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="AY293618.1_Gavia_stellata_mitochondrion"]="Neoaves"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="FJ769841.1_Balearica_regulorum_mitochondrion"]="Neoaves"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="AP010797.1_Accipiter_gentilis_mitochondrial_DNA"]="Neoaves"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="MT920475.1_Cacatua_alba_mitochondrion"]="Neoaves"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="MF431746.1_Strix_occidentalis_caurina_voucher"]="Neoaves"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="KP319029.1_Columba_livia_mitochondrion"]="Neoaves"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="CM016612.1_Calypte_anna_isolate_BGI_N300_mitochondrion"]="Neoaves"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="NC_045372.1_Picus_canus_mitochondrion"]="Neoaves"

MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="AP003322.1_Gallus_gallus_gallus_mitochondrial_DNA"]="Galloanserea"
MergedNUMTtable$Clade[MergedNUMTtable$nMITO=="AF338715.1_Struthio_camelus_mitochondrion"]="Palaeognathae"

MergedNUMTtable$SuperClade=NA
MergedNUMTtable$SuperClade[MergedNUMTtable$Clade=="Kestrel"]="Falco"
MergedNUMTtable$SuperClade[MergedNUMTtable$Clade=="Mid-Falcon"]="Falco"
MergedNUMTtable$SuperClade[MergedNUMTtable$Clade=="Hierofalco"]="Falco"
MergedNUMTtable$SuperClade[MergedNUMTtable$Clade=="Hierofalco"]="Falco"
MergedNUMTtable$SuperClade[MergedNUMTtable$Clade=="Peregrine"]="Falco"
MergedNUMTtable$SuperClade[MergedNUMTtable$Clade=="Mid-Falcon"]="Falco"
MergedNUMTtable$SuperClade[MergedNUMTtable$Clade=="sparverius"]="Falco"
MergedNUMTtable$SuperClade[MergedNUMTtable$Clade=="Polyborinae"]="Falconidae"
MergedNUMTtable$SuperClade[MergedNUMTtable$Clade=="Herpetherinae"]="Falconidae"
MergedNUMTtable$SuperClade[MergedNUMTtable$Clade=="Neoaves"]="Neoaves"
MergedNUMTtable$SuperClade[MergedNUMTtable$Clade=="Galloanserea"]="Galloanserea"
MergedNUMTtable$SuperClade[MergedNUMTtable$Clade=="Palaeognathae"]="Palaeognathae"

MergedNUMTtable$CladeCount=NA
MergedNUMTtable.numt=MergedNUMTtable[MergedNUMTtable$TYPE=="NUMT",]
MergedNUMTtable.numt.purge=MergedNUMTtable.numt[which(!is.na(MergedNUMTtable.numt$BASEPAIRS)),]
for (c in unique(MergedNUMTtable.numt.purge$Clade)){
MergedNUMTtable.numt.purge$CladeCount[MergedNUMTtable.numt.purge$Clade==c]=length(which(MergedNUMTtable.numt.purge$Clade==c))
}

setwd("~/NYU_AD/FalconProject/Notebooks/FalconSeq_and_Reseq2018/RefGenomeAnalysis/Synthesis/Data/ChromAnnotated")
Isochores100KB=read.delim("AllultChromAnnotated100KBWindowsSummary.txt", header=TRUE, sep="\t")
Isochores100KB$WinPos=trunc(as.numeric(gsub(".*__([0-9]+)__.*","\\1",Isochores100KB$Window))/100000)
MergedNUMTtable.numt.purge$Genome=gsub(":.*","", row.names(MergedNUMTtable.numt.purge))
MergedNUMTtable.numt.purge$genScaff=gsub("_.*","_1",gsub("[a-zA-Z]+[0-9]*:","", row.names(MergedNUMTtable.numt.purge)))
MergedNUMTtable.numt.purge$WindowNumStart=trunc(MergedNUMTtable.numt.purge$genPOS1/100000)
MergedNUMTtable.numt.purge$WindowNumEnd=trunc(MergedNUMTtable.numt.purge$genPOS2/100000)
MergedNUMTtable.numt.purge$ancChromWindowStart="NA"
MergedNUMTtable.numt.purge$ancChromWindowEnd="NA"
MergedNUMTtable.numt.purge$FalcoChromWindowStart="NA"
MergedNUMTtable.numt.purge$FalcoChromWindowEnd="NA"
for (n in row.names(MergedNUMTtable.numt.purge[MergedNUMTtable.numt.purge$WindowNumStart>=1,])){
Temp=Isochores100KB[Isochores100KB$Genome==MergedNUMTtable.numt.purge[n,]$Genome,]
Temp=Temp[Temp$REF_SCAFFOLD==MergedNUMTtable.numt.purge[n,]$genScaff,]
if (MergedNUMTtable.numt.purge[n,]$WindowNumStart<=max(Temp$WinPos) & MergedNUMTtable.numt.purge[n,]$WindowNumEnd<=max(Temp$WinPos)){
MergedNUMTtable.numt.purge[n,]$ancChromWindowStart=as.character(Temp[Temp$WinPos==MergedNUMTtable.numt.purge[n,]$WindowNumStart,]$CONSENSUS_CLASS_MATCH)
MergedNUMTtable.numt.purge[n,]$ancChromWindowEnd=as.character(Temp[Temp$WinPos==MergedNUMTtable.numt.purge[n,]$WindowNumEnd,]$CONSENSUS_CLASS_MATCH)
MergedNUMTtable.numt.purge[n,]$FalcoChromWindowStart=as.character(Temp[Temp$WinPos==MergedNUMTtable.numt.purge[n,]$WindowNumStart,]$FalcoCONSENSUS_CLASS_MATCH)
MergedNUMTtable.numt.purge[n,]$FalcoChromWindowEnd=as.character(Temp[Temp$WinPos==MergedNUMTtable.numt.purge[n,]$WindowNumStart,]$FalcoCONSENSUS_CLASS_MATCH)
}
}

MergedNUMTtable.numt.purge.filt=MergedNUMTtable.numt.purge[MergedNUMTtable.numt.purge$ancChromWindowStart==MergedNUMTtable.numt.purge$ancChromWindowEnd & MergedNUMTtable.numt.purge$FalcoChromWindowStart==MergedNUMTtable.numt.purge$ancChromWindowEnd,]
MergedNUMTtable.numt.purge.filt.noamb=MergedNUMTtable.numt.purge[MergedNUMTtable.numt.purge$ancChromWindowStart!="Ambiguous" & MergedNUMTtable.numt.purge$FalcoChromWindowStart!="Ambiguous",]
MergedNUMTtable.numt.purge.filt.noamb.falco=MergedNUMTtable.numt.purge.filt.noamb[MergedNUMTtable.numt.purge.filt.noamb$nSuperClade=="Falco",]
FalcoChromLarge=dim(MergedNUMTtable.numt.purge.filt.noamb.falco[MergedNUMTtable.numt.purge.filt.noamb.falco$FalcoChromWindowStart!="Microchromosome",])[1]
FalcoChromMicro=dim(MergedNUMTtable.numt.purge.filt.noamb.falco[MergedNUMTtable.numt.purge.filt.noamb.falco$ChromWindowStart=="Microchromosome",])[1]

AncChromLarge=dim(MergedNUMTtable.numt.purge.filt.noamb[MergedNUMTtable.numt.purge.filt.noamb$ancChromWindowStart!="Microchromosome",])[1]
AncChromMicro=dim(MergedNUMTtable.numt.purge.filt.noamb[MergedNUMTtable.numt.purge.filt.noamb$ancChromWindowStart=="Microchromosome",])[1]

FalcoChromLarge=dim(MergedNUMTtable.numt.purge.filt.noamb[MergedNUMTtable.numt.purge.filt.noamb$SuperClade=="Falco" & MergedNUMTtable.numt.purge.filt.noamb$FalcoChromWindowStart!="Microchromosome",])[1]
FalcoChromMicro=dim(MergedNUMTtable.numt.purge.filt.noamb[MergedNUMTtable.numt.purge.filt.noamb$SuperClade=="Falco" & MergedNUMTtable.numt.purge.filt.noamb$FalcoChromWindowStart=="Microchromosome",])[1]
PropFalcoMicroNUMT=FalcoChromMicro/FalcoChromLarge
Isochores100KB.Falco_noamb=Isochores100KB[Isochores100KB$FalcoCONSENSUS_CLASS_MATCH!="Ambiguous",]
FalcoIsochoresLarge=dim(Isochores100KB.Falco_noamb[Isochores100KB.Falco_noamb$FalcoCONSENSUS_CLASS_MATCH!="Microchromosome",])[1]
FalcoIsochoresMicro=dim(Isochores100KB.Falco_noamb[Isochores100KB.Falco_noamb$FalcoCONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
PropFalcoIsochroesMicro=FalcoIsochoresMicro/FalcoIsochoresLarge
binom.test(FalcoChromMicro,FalcoChromLarge+FalcoChromMicro,p=PropFalcoIsochroesMicro, alternative="two.sided")
#proportion NUMTs inserting on Microchromosomes in Falcons followed by proportion of Isochores on microchromosomes
PropFalcoMicroNUMT
PropFalcoIsochroesMicro

AncChromLarge=dim(MergedNUMTtable.numt.purge.filt.noamb[MergedNUMTtable.numt.purge.filt.noamb$SuperClade!="Falco" & MergedNUMTtable.numt.purge.filt.noamb$FalcoChromWindowStart!="Microchromosome",])[1]
AncChromMicro=dim(MergedNUMTtable.numt.purge.filt.noamb[MergedNUMTtable.numt.purge.filt.noamb$SuperClade!="Falco" & MergedNUMTtable.numt.purge.filt.noamb$FalcoChromWindowStart=="Microchromosome",])[1]
PropAncMicroNUMT=AncChromMicro/AncChromLarge
Isochores100KB.anc_noamb=Isochores100KB[Isochores100KB$CONSENSUS_CLASS_MATCH!="Ambiguous",]
AncIsochoresLarge=dim(Isochores100KB.anc_noamb[Isochores100KB.Falco_noamb$CONSENSUS_CLASS_MATCH!="Microchromosome",])[1]
AncIsochoresMicro=dim(Isochores100KB.anc_noamb[Isochores100KB.Falco_noamb$CONSENSUS_CLASS_MATCH=="Microchromosome",])[1]
PropAncIsochroesMicro=AncIsochoresMicro/AncIsochoresLarge
binom.test(AncChromMicro,AncChromLarge+AncChromMicro,p=PropAncIsochroesMicro, alternative="two.sided")
#proportion NUMTs inserting on Microchromosomes in Ancesters of Falcons followed by proportion of Isochores on microchromosomes
PropAncMicroNUMT
PropAncIsochroesMicro

fisher.test(as.matrix(rbind(cbind(FalcoChromMicro, FalcoChromLarge), cbind(AncChromMicro, AncChromLarge))))
binom.test(FalcoChromMicro,FalcoChromLarge+FalcoChromMicro,p=PropAncMicroNUMT, alternative="two.sided")

MergedNUMTtable$CladeCount=NA
MergedNUMTtable.numt=MergedNUMTtable[MergedNUMTtable$TYPE=="NUMT",]
MergedNUMTtable.numt.purge=MergedNUMTtable.numt[which(!is.na(MergedNUMTtable.numt$BASEPAIRS)),]
for (c in unique(MergedNUMTtable.numt.purge$Clade)){
MergedNUMTtable.numt.purge$CladeCount[MergedNUMTtable.numt.purge$Clade==c]=length(which(MergedNUMTtable.numt.purge$Clade==c))
}
print(unique(paste(MergedNUMTtable.numt.purge$Clade, MergedNUMTtable.numt.purge$CladeCount)))

summary.lm(glm(MergedNUMTtable.numt.purge$BASEPAIRS~MergedNUMTtable.numt.purge$sMitoDist, family="poisson"))

cladeBP.aov=aov(log2(MergedNUMTtable.numt.purge$BASEPAIRS)~as.factor(MergedNUMTtable.numt.purge$Clade))
cladeBP.aov.tukey=TukeyHSD(cladeBP.aov, conf.level=0.95)
summary.lm(cladeBP.aov)
cladeBP.aov.tukey

summary.lm(glm(MergedNUMTtable.numt.purge$BASEPAIRS~MergedNUMTtable.numt.purge$sMitoDist, family="poisson"))
cladeBP.aov=glm(log2(MergedNUMTtable.numt.purge$BASEPAIRS)~as.factor(MergedNUMTtable.numt.purge$Clade), family="poisson")
summary.lm(cladeBP.aov)
scladeBP.aov=glm(MergedNUMTtable.numt.purge$BASEPAIRS~as.factor(MergedNUMTtable.numt.purge$SuperClade), family="poisson")
summary.lm(scladeBP.aov)
pairwise.wilcox.test(MergedNUMTtable.numt.purge$BASEPAIRS, as.factor(MergedNUMTtable.numt.purge$SuperClade), p.adjust.method ="bonferroni")
aggregate(MergedNUMTtable.numt.purge$BASEPAIRS, list(MergedNUMTtable.numt.purge$SuperClade), FUN=mean)

pairwise.wilcox.test(MergedNUMTtable.numt.purge$BASEPAIRS, as.factor(MergedNUMTtable.numt.purge$Clade), p.adjust.method ="bonferroni")
aggregate(MergedNUMTtable.numt.purge$BASEPAIRS, list(MergedNUMTtable.numt.purge$Clade), FUN=mean)

library(ggplot2)
MergedNUMTtable.numt.purge$Clade=factor(MergedNUMTtable.numt.purge$Clade, levels=c("Peregrine","Hierofalco", "Mid-Falcon", "Kestrel", "Polyborinae","Herpetherinae","Neoaves", "Galloanserea", "Palaeognathae" ))
ggplot(MergedNUMTtable.numt.purge, aes(x=MergedNUMTtable.numt.purge$Clade, y=MergedNUMTtable.numt.purge$BASEPAIRS))+geom_boxplot()+xlab("Clade")+labs(y='NUMT Length (BP)')+theme(axis.text.y=element_text(size=7),axis.text.x=element_text(size=7), axis.title=element_text(size=12, face="plain"))+theme(strip.text = element_text(face = "bold"))
ggsave("Figure5_NUMTsLength-by-Node.png", plot = last_plot(), device = png(), path = "~/NYU_AD/FalconProject/Manuscripts/FalconRef8/Figures/Ultimate/Intermediate", scale = 1, width = 7, height = 5.0, units = c("in"), dpi =1200, limitsize = TRUE)

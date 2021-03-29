#Assignment 2 Project 4

# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

#set your working directory
setwd("~/Desktop/BU/PhD/Spring_2021_classes/Ecological_genomics/Assign2/Project4") #you will need to change to your own directory

#Bioconductor install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

###conduct array quality metrics to detect and remove outliers
library(DESeq2) 
library(affycoretools)
library(arrayQualityMetrics)
library(genefilter)
library(Biobase)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(vegan)
library(ggrepel)
library(tidyverse)

#read in iso2gene tab separated file
gg=read.table("Bmin_iso2gene.tab",sep="\t")

#read in counts 
countData <- read.table("allcounts_Harvey_Sym.txt") #what do these 3 groups represent?
head(countData)
length(countData[,1])

#Renaming columns, which I don't think we need to do??
#names(countData)=c( "pH7.5a", "pH7.5b", "pH7.5c", "pH7.6a", "pH7.6b", "pH7.6c", "pH8a", "pH8b",  "pH8c")
#row.names(countData)=sub("", "isogroup", rownames(countData))
#head(countData)


# # # #look for outliers
treat=c( "Recovery", "Recovery", "Recovery", "Stress", "Stress", "Stress") #creating a condition table
g=data.frame(treat)
g
colData= g

#ASK WHAT THIS IS
#converting counts data, treatment info, and statistical model we want to use (design) into a DESeq object so we can work w it in DESeq
dds=DESeqDataSetFromMatrix(countData=countData,
                           colData = g,
                           design = ~treat)

vsd.ge=assay(vst(dds))
rl=vst(dds)
e=ExpressionSet(assay(rl), AnnotatedDataFrame(as.data.frame(colData(rl))))
arrayQualityMetrics(e, intgroup=c("treat"),force=T) #will tell us about outliers in our data

# double-click index.html
 #go to 1:19 on class recording to see explanation of the index.html

 #####you only ever need to run the above code once. Outliers are decided at the beginning. 
## I like to close R and restart with packages etc
##So, please save your script and restart R - Sarah

#fav_stressC detected as outliers via the Distances between arrays method (Figure 2) and MA plot method (Figure 9). For the purpose of this assignment, we will include it, but in a real scenario we would exclude it.

library("DESeq2")
library("ggplot2")

#one step DESeq - only 9 samples here, with large datasets could take a very long time, may want to use SCC
dds<-DESeq(dds)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

head(dds)
res<- results(dds) #these are our results
res
#just do this to make sure there's nothing insanely funky going on 
#post-normalized counts - pvalue corresponds to see how differentially expressed a gene is between treatments
#padj adjusts for mutliple tests
#baseMean is avg count for one comparison (whichever treatment group came first)
#log2FoldChange - diff in counts btwn basemean vs other treatment group
#symcomp10000 (this is a gene) - the recovery treatment group had avg of 5.57 counts; stress group had a 2fold change (logtransformed) 

#Look at dispersions plot
plotDispEsts(dds, main="Wright et al. 2019 symbiont dispersion plot")
#should look like hockey stick
#this is the normalization method
#visual of what DEseq2 is doing to your data; fitting it to this curve

####################Stress vs recovery pairwise comparisons
#can only do pairwise comparisons
colData$Stress<-factor(colData$treat, levels=c("Stress","Recovery")) #levels matter bc gene expression is always relative; second term is the one that's relative (treatment 1st, control 2nd)
##second term is the "control"
resStress <- results(dds, contrast=c("treat","Stress","Recovery"))
#how many FDR < 10%?
table(resStress$padj<0.01) #5093 reads, p-adjusted reads are 6 --> always use p-adjusted to look at numbers of differentially expressed genes - here we have 6 diferentially expressed genes
#p-adjusted accounts for all the different times you're comparing things - could get false positive by chance without adjustment
# 0.1=486
# 0.05=381
# 0.01=242
summary(resStress) #lots of genes removed bc of low counts; 0.1% of genes upregulated in stress relative to recovery, and 0.43% were downregulated under stress relative to recovery

#another way to look at it, no. of differentially expressed genes (6)
nrow(resStress[resStress$padj<0.01 & !is.na(resStress$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   #228

plotMA(resStress, main="Stress vs Recovery")
plotMA(resStress, main="Stress vs Recovery", ylim=c(-2,2)) #sometimes nice to have a ylim if your axes are weird
#MA plot shows mean of normalized counts; left is more lowly expressed genes, right is higher expressed genes
#fewer differentially expressed genes (in red) upregulated ; most of genes are downregulated (blue?)
#ASK ABOUT THIS -- what do the colors mean again?? is blue downregulated? red upregulated and we just have so few we can't see them? grey not differentially expressed?

#put results into a dataframe
results <- as.data.frame(resStress)
head(results)

nrow(resStress[resStress$padj<0.1 & resStress$log2FoldChange > 0 & !is.na(resStress$padj),]) #no. genes upregulated
nrow(resStress[resStress$padj<0.1 & resStress$log2FoldChange < 0 & !is.na(resStress$padj),]) #no. genes downregulated
#23 upregulated
#98 downregulated

#sometimes include this table as supp info in a paper
write.table(resStress, file="StressvRecovery.txt", quote=F, sep="\t")

cd <- read.table("StressvRecovery.txt")
head(cd)

##make the GO table for MWU
head(cd)

library(dplyr)
cd
go_input = cd %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
#above code just took negative log of p value, then making either positive (upregulated) or negative (downregulated); ranked p value that also has directionality
head(go_input) #one row for each gene that was diferentially expressed
colnames(go_input) <- c("gene", "pval")
head(go_input)
write.csv(go_input, file="Stress_GO.csv", quote=F, row.names=FALSE) #this will be input for GO


###############################################################################################
##############################################################################
#--------------get pvals
valStress=cbind(resStress$pvalue, resStress$padj)
head(valStress)
colnames(valStress)=c("pval.Stress", "padj.Stress")
length(valStress[,1])
table(complete.cases(valStress))

valStress=cbind(resStress$pvalue, resStress$padj)
head(valStress)
colnames(valStress)=c("pval.Stress", "padj.Stress")
length(valStress[,1])
table(complete.cases(valStress))

######-------------make rlogdata and pvals table
rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
colnames(rld)=paste(colData$treat)
head(rld)
length(rld[,1])

rldpvals=cbind(rld,valStress, valStress)
head(rldpvals)
dim(rldpvals)
# [1] 19717    13
table(complete.cases(rldpvals))
#FALSE  TRUE 
#17202  2515 

write.csv(rldpvals, "RLDandPVALS.csv", quote=F)

colnames(rld)=paste(colData$treat)
head(rld)

library(RColorBrewer)
# Sample distance heatmap
sampleDists <- as.matrix(dist(t(rld)))
library(gplots)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          margin=c(10, 10), main="Sample Distance Matrix")

# 2 stress samples that act similarly, and one stress sample that has different gene expression...lots of genes doing same thing, or different samples respond differently
#if there's a strong signature, all stress samples would cluster together..here, 2/3 stress samples cluster, but third is odd


#################################################################################
###########################heat map of sample distances for pco2
rldpvals <- read.csv(file="Crep2016_RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld=rldpvals[,1:9]
head(rld)

sampleDists <- dist(t(rld))
sampleDistMatrix <- as.matrix( sampleDists )
treat=c( "pH7.5", "pH7.5", "pH7.5", "pH7.6", "pH7.6", "pH7.6", "pH8", "pH8",  "pH8")
colnames(sampleDistMatrix)=paste(treat)
rownames(sampleDistMatrix)=paste(treat)

library("pheatmap")
heat.colors = colorRampPalette(rev(c("blue","yellow","red")),bias=0.3)(100)
pheatmap(sampleDistMatrix,color = heat.colors,cex=0.9,border_color=NA,cluster_rows=T,cluster_cols=T)

library(vegan)
library(ggplot2)
library(ggrepel)
library(tidyverse)

rld_t=t(rld)
pca <- prcomp(rld_t,center = TRUE, scale. = TRUE)
head(pca)
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = row.names(pca_s)
pca_s$treat=colData$treat
head(pca_s)

cbPalette <- c("darkgoldenrod2",  "darkolivegreen3", "dodgerblue3")
ggplot(pca_s, aes(PC1, PC2, color = treat, pch = treat)) +
  geom_point(size=3) +
  #  geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  # geom_density2d(alpha=.5)+
  geom_polygon(alpha=.2)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
  head(pca)
adonis(pca$x ~ treat, data = pca_s, method='eu', na.rm = TRUE)
          # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treat      2     40140 20069.9  13.048 0.81306  0.003 **
# Residuals  6      9229  1538.1         0.18694          
# Total      8     49369                 1.00000 


###################################heatmaps for genes NS vs FR
rldpvals <- read.csv(file="Crep2016_RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld_site= rldpvals[,1:9]
head(rld_site)
gg=read.table("Crep454_iso2gene.tab",sep="\t", row.names=1)
head(gg)

nrow(rldpvals[rldpvals$padj.Stress<0.01& !is.na(rldpvals$padj.76),])
#242

topnum= 100 # number of DEGS
head(rldpvals)
top100=head(rldpvals[order(rldpvals$padj.76), ],topnum)
head(top100)
length(top100[,1])
summary(top100)
###
library(pheatmap)
head(top100)
p.val=0.1 # FDR cutoff
conds=top100[top100$padj.76<=p.val & !is.na(top100$padj.76),]
length(conds[,1])

exp=conds[,1:9] # change numbers to be your vsd data columns
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them
head(explc)

ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

pheatmap(explc,cluster_cols=T,scale="row",color=col0, show_rownames = F)

###################################Heatmap for the genes in common
rldpvals <- read.csv(file="Crep2016_RLDandPVALS.csv", row.names=1)
head(rldpvals)
p.val=0.1 # FDR cutoff
conds=rldpvals[rldpvals$padj.76<=p.val & !is.na(rldpvals$padj.76) & rldpvals$padj.Stress<=p.val & !is.na(rldpvals$padj.Stress),]
rld_data= conds[,c(1:9)]
head(rld_data)
nrow(rld_data)
gg=read.table("Crep454_iso2gene.tab",sep="\t", row.names=1)
library(pheatmap)
means=apply(rld_data,1,mean) # means of rows
explc=rld_data-means # subtracting them

ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

pheatmap(explc,cluster_cols=T,scale="row",color=col0, show_rownames = F)

# Make annotation table for pheatmap
ann = data.frame(cond = c('7.5', '7.5', '7.5', '7.6', '7.6', '7.6', '8', '8', '8'))
rownames(ann) <- names(explc)

# Set colors
Var1        <- c("darkgoldenrod2",  "darkolivegreen3", "dodgerblue3")
names(Var1) <- c("7.5", "7.6", "8")
anno_colors <- list(cond = Var1)

pheatmap(as.matrix(explc),annotation_col=ann,annotation_colors=anno_colors,cex=1.2,color=col0,border_color=NA,clustering_distance_rows="correlation",clustering_distance_cols="correlation", show_rownames=T)

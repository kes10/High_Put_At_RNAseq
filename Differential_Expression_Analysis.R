### Appendix E: Differential Expression Analysis 

## This goes through the generation of a count matrix with Feature Counts
## Then DE analysis with edgeR including filtering, normalization,
## DMS plots, treatment contrasts, Top tag lists, and a heatmap

#Install necessary packages


#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#biocLite('GenomicRanges')
#biocLite('GenomicFeatures')
#biocLite('Rsamtools')
#biocLite("org.Mm.eg.db")
#source("http://bioconductor.org/biocLite.R")
#biocLite("Rsubread")

library('edgeR')
library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(Rsubread)

## moved my Bam files (post mapping) to a usable location on my own desktop. It was
## at this stage that I brought together the data from both runs by merging and sorting
## the bam files.

## I used the featureCounts program from Rsubread to generate a count matrix (as recommended
## in the edgeR manual). It needs annotation gtf file (have from the Ensembl build,
## and the bam files for the libraries after mapping)

setwd("~/Desktop/RNA_seq/Analysis/featureCounts_BothLanes")
fc<-featureCounts(annot.ext="genes.gtf", isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_id", files=c("Cnst0_1_fc.srt.bam","Cnst0_2_fc.srt.bam","Cnst0_3_fc.srt.bam","Cnst24_1_fc.srt.bam","Cnst24_2_fc.srt.bam","Cnst24_3_fc.srt.bam","Ind0_1_fc.srt.bam","Ind0_2_fc.srt.bam","Ind0_3_fc.srt.bam","Indi24_1_fc.srt.bam","Indi24_2_fc.srt.bam","Indi24_3_fc.srt.bam","Indi48_1_fc.srt.bam","Indi48_2_fc.srt.bam","Indi48_3_fc.srt.bam","Indn24_1_fc.srt.bam","Indn24_2_fc.srt.bam","Indn24_3_fc.srt.bam","Indn48_1_fc.srt.bam","Indn48_2_fc.srt.bam","Indn48_3_fc.srt.bam","Wt0_1_fc.srt.bam","Wt0_2_fc.srt.bam","Wt0_3_fc.srt.bam"), isPairedEnd=TRUE, useMetaFeatures=TRUE)

#PUT COUNTS INTO TABLE
write.table(x=data.frame(fc$annotation[,c("GeneID","Length")],fc$counts,stringsAsFactors=FALSE),file="counts_BL.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(x=data.frame(fc$annotation[,c("GeneID")],fc$counts,stringsAsFactors=FALSE),file="counts_BL2.txt",quote=FALSE,sep="\t",row.names=FALSE)
counts_BL <- read.delim("~/Desktop/RNA_seq/Analysis/featureCounts_BothLanes/counts_BL.txt", header=TRUE)
head(counts_BL)
count.dataBL <-counts_BL[,-c(1,ncol(counts))]    #not sure what this does, I think removes the gene as a column so you can add it in as a dimension name. 
head(count.dataBL)
rownames(count.dataBL) <-counts_BL[,1] #gene names as row dim names
head(count.dataBL)
colnames(count.dataBL) <-paste(c(rep("Length",1),rep("Cnst0",3),rep("Cst24",3),rep("Indl0",3),rep("Idi24",3),rep("Idi48",3),rep("Idn24",3),rep("Idn48",3),rep("Wldt0",3)),c(1:1,1:3,1:3,1:3,1:3,1:3,1:3,1:3,1:3),sep=" ") #sample names
head(count.dataBL)

#SOME BASIC INFO
summary(count.dataBL)
dim(count.dataBL)    #33602 and 25 rows
colSums(count.dataBL)     #library sizes
colSums(count.dataBL)/1e06 #library sizes in millions of reads
table(rowSums(count.dataBL))[1:30]    #number of genes with low counts (counts 36 -> 92)

##CREATE GROUP SET UP LIST. Will give this to DGEList along with the count matrix and then make the edgeR object for analysis
count.dataBL2 <- count.dataBL[,-1]  #get rid of the length information for the DE analysis.
head(count.dataBL2)

#SET UP EXPERIMENTAL DESIGN
Treat <- factor(substring(colnames(count.dataBL2),1,5)) #takes the column names first 4 spaces/letters
Treat <- relevel(Treat, ref="Wldt0") #sets Wild type as my reference (control) and places it first in the treatment line up (8 levels: Wlt0 Cst0 Cst2 Ii24 Ii48 In24 In48 Ind0)
Block <- factor(substring(colnames(count.dataBL2),7,7))  #takes space 6 of column names only
Treat 
Block

#FILTER DATA 
dim(count.dataBL2)   # 33602 genes for 24 samples
keep <- rowSums(cpm(count.dataBL2)>5) >=3    #filtering - requires a gene to have at least 5 cpm in at least 3 of the samples
A <- count.dataBL2[keep, ]
dim(A)  #for filter of at least 2 cpm for 3 samples = 19068, for filter of 5 = 16862

#CREATE DGEList and apply TMM normalization
#DATA NORMALIZATION: edgeR normalizes by total count. finds set of scaling factors for lib sizes that minimize log-fold changes between samples for most genes. default method is to use trimmed mean of M-values (TMM)
DGE <- DGEList(counts=A,group=Treat)
DGE$samples  #for instance, Cnst0_1 = 15965960, all norm factors = 1
DGE <- calcNormFactors(DGE)
DGE$samples$lib.size<-colSums(DGE$counts) # resets the lib size now that some genes have been lost 
DGE$samples #now Cnst0 norm factor = 1.03 etc

#SET UP MODEL -DESIGN MATRIX
designA <- model.matrix(~Block+Treat)
rownames(designA) <- colnames(DGE)
designA
designB <- model.matrix(~0+Block+Treat) #including the zero in the matrix allows a space for the control (Wt0). doesn't change anything, just the visual
rownames(designB) <- colnames(DGE)
colnames(designA)
designB
designC <-model.matrix(~0+Treat) #not including the blocking into design
designC

# DISPERSION ESTIMATES
DGE <- estimateGLMCommonDisp(DGE, designA, verbose = T)  #average over all genes, output: Disp = 0.06888, BCV = 0.2624 with Design A 
DGE <- estimateGLMTrendedDisp(DGE, designA) # allows for a possible abudance trend
DGE <- estimateGLMTagwiseDisp(DGE, designA) #per tag (per gene?)
plotBCV(DGE)

#MDS PLOT: look at how samples relate based on multidimensional scaling. Samples closer to each other are more similar and vice versa.
plotMDS(DGE, method="bcv", col=as.numeric(DGE$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:8, pch=20)
  # -> Definitely still see some effects of time and also some block effects 
  
  
###TREATMENT PAIRWISE COMPARISONS, just changed the index in the lrt <- statement and went through them manually.
alltreat.pair.contrast.matrix = cbind(C0C24 = c(0,0,0,1,-1,0,0,0,0), C0Ii24 = c(0,0,0,1,0,-1,0,0,0,0),
                                      C24Ii24 = c(0,0,0,0,1,-1,0,0,0,0), C0Ii48 = c(0,0,0,1,0,0,-1,0,0,0),
                                      C24Ii48 = c(0,0,0,0,1,0,-1,0,0,0), Ii24Ii48 = c(0,0,0,0,0,1,-1,0,0,0),
                                      C0In24 = c(0,0,0,1,0,0,0,-1,0,0), C24In24 = c(0,0,0,0,1,0,0,-1,0,0),
                                      Ii24In24 = c(0,0,0,0,0,1,0,-1,0,0), Ii48In24 = c(0,0,0,0,0,0,1,-1,0,0),
                                      C0In48 = c(0,0,0,1,0,0,0,0,-1,0), C24In48 = c(0,0,0,0,1,0,0,0,-1,0),
                                      Ii24In48 = c(0,0,0,0,0,1,0,0,-1,0), Ii48In48 = c(0,0,0,0,0,0,1,0,-1,0),
                                      In24In48 = c(0,0,0,0,0,0,0,1,-1,0), C0I0 = c(0,0,0,1,0,0,0,0,0,-1),
                                      C24I0= c(0,0,0,0,1,0,0,0,0,-1), Ii24I0 = c(0,0,0,0,0,1,0,0,0,-1),
                                      Ii48I0 = c(0,0,0,0,0,0,1,0,0,-1), In24I0 = c(0,0,0,0,0,0,0,1,0,-1),
                                      In48I0 = c(0,0,0,0,0,0,0,0,1,-1))
  
fitA <- glmFit(DGE, designA)
fitB <- glmFit(DGE, designB)
fitC <- glmFit(DGE, designC)

lrt <- glmLRT(fitB, contrast = alltreat.pair.contrast.matrix[,1])  #shift through all the contrast statements in the contrast matrix
FDR <- p.adjust(lrt$table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrt))  
topTags(lrt)
#for wild type comparisons (because WT is my reference level) :
lrt <- glmLRT(fitA, coef = 10)
FDR <- p.adjust(lrt$table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrt))  
topTags(lrt)

#collected up and down-regulated DE gene #s in a csv file for plotting
#Import data for ALL PAIRWISE COMPARISONS GRAPH
allpairwise_DE2s <- read.csv("~/Desktop/RNA_seq/Analysis/featureCounts_BothLanes/allpwc_bothlanes.csv")  #version to has no comparison names so the lableing doesn't get messed up
View(allpairwise_DE2s)
head(allpairwise_DE2s)
attach(allpairwise_DE2s)
levels(Comparison)
levels(Direction)
data = DE_Genes
head(data)
data = matrix(data, ncol = 28, byrow = F)
colnames(data) = levels(Comparison)
rownames(data) = levels(Direction)
par(mar = c(5.1, 4.1, 4.1, 7.1), xpd = T) #play around with the margins so that adding in the legend doesn't mess everything up
rgb.palette1 <- colorRampPalette(c("white", "black"), space = "Lab")  #switch around coloring if you like

#Make a percentage based stacked bar plot (make sure ordering is right here)
prop = prop.table(data, margin = 2)
barplot(prop, col = heat.colors(length(rownames(prop))), width = 2)
legend("center", inset = c(-0.05,0.05), fill = heat.colors(length(rownames(prop))), legend = rownames(data))

#for counts based vertical stacked barplot
par(mar = c(5.1, 4.1, 5.1, 2.1), xpd = T)
barplot(data,col = rgb.palette1(length(rownames(data))), space = 0.5, ylab = "DE Genes", ylim= c(0,20000), xlab = "   ")
legend("center", inset = c(-0.000005,0.05), fill = rgb.palette1(length(rownames(prop))), legend = c("Up","Down","Unchanged"))
  ## change axis labels to be slanted
x <-barplot(data,col = rgb.palette1(length(rownames(data))), space = 0.5, ylab = "DE Genes", ylim= c(0,20000), xlab = "")
text(x, par("usr")[3]-1, srt = 60, adj = 1.2, labels = c("C0C24","C0Ii24","C24Ii24","C0Ii48","C24Ii48","Ii24Ii48","C0In24","C24In24","Ii24In24","Ii48In24","C0In48","C24In48","Ii24In48","Ii48In48","In24In48","C0I0","C24I0","Ii24I0","Ii48I0","In24I0","In48I0","C0W0","C24W0","Ii24W0","Ii48W0","In24W0","In48W0","I0W0"),xpd = TRUE, font = 0.1)

## MAJOR CONTRAST Q's
lrtCvsIi1 <- glmLRT(fitA, contrast = c(0,0,0,1,1,-1,-1,0,0,0))
FDR <- p.adjust(lrtCvsIi1$table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrtCvsIi1))  
results1 <- topTags(lrtCvsIi1, n=893)
head(results1)
write.table(as.matrix(results$table),file="DE_CvsIi1.txt",sep="\t")
dt.CvsIitags1 <- rownames(DGE[as.logical(dt)])
plotSmear(lrtCvsIi1, de.tags=dt.CvsIitags1)
abline(h=c(-2,2), col = "purple")

lrtCvsW <- glmLRT(fitA, contrast = c(0,0,0,-0.5,-0.5,0,0,0,0,0))  #not 100% sure this is right...
FDR <- p.adjust(lrtCvsW$table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrtCvsW))  
results <- topTags(lrtCvsW, n=293)
head(results)
write.table(as.matrix(results$table),file="DE_CvsW.txt",sep="\t")
dt.CvsWtags <- rownames(DGE[as.logical(dt)])
plotSmear(lrtCvsW, de.tags=dt.CvsWtags)
abline(h=c(-2,2), col = "purple")

lrtIivsIn <- glmLRT(fitA, contrast = c(0,0,0,0,0,1,1,-1,-1,0)) 
FDR <- p.adjust(lrtIivsIn$table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrtIivsIn))  
topTags(lrtIivsIn)
results <- topTags(lrtIivsIn, n=1)
write.table(as.matrix(results$table),file="DE_IivsIn.txt",sep="\t")
dt.IivsIntags <- rownames(DGE[as.logical(dt)])
plotSmear(lrtIivsIn, de.tags=dt.IivsIntags)
abline(h=c(-2,2), col = "purple")
lrtIivsIn0 <- glmLRT(fitA, contrast = c(0,0,0,0,0,3,3,-2,-2,-2)) 
FDR <- p.adjust(lrtIivsIn0$table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrtIivsIn0))  
topTags(lrtIivsIn0)
results <- topTags(lrtIivsIn0, n=1)
write.table(as.matrix(results$table),file="DE_IivsIn0.txt",sep="\t")
dt.IivsIn0tags <- rownames(DGE[as.logical(dt)])
plotSmear(lrtIivsIn0, de.tags=dt.IivsIn0tags)
abline(h=c(-2,2), col = "purple")

lrtWvsIn <- glmLRT(fitA, contrast = c(0,0,0,0,0,0,0,-(1/3),-(1/3),-(1/3))) 
FDR <- p.adjust(lrtWvsIn$table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrtWvsIn))  
topTags(lrtWvsIn)
results <- topTags(lrtWvsIn, n=1202)
write.table(as.matrix(results$table),file="DE_WvsIn.txt",sep="\t")
dt.WvsIntags <- rownames(DGE[as.logical(dt)])
plotSmear(lrtWvsIn, de.tags=dt.WvsIntags)
abline(h=c(-2,2), col = "purple")

lrt0vs24 <- glmLRT(fitA, contrast = c(0,0,0,1,-1,0,0,-1,0,1)) 
FDR <- p.adjust(lrt0vs24$table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrt0vs24))  
topTags(lrt0vs24)
results <- topTags(lrt0vs24, n=882)
write.table(as.matrix(results$table),file="DE_0vs24.txt",sep="\t")
dt.0vs24tags <- rownames(DGE[as.logical(dt)])
plotSmear(lrt0vs24, de.tags=dt.0vs24tags)
abline(h=c(-2,2), col = "purple")

lrt24vs48 <- glmLRT(fitA, contrast = c(0,0,0,0,0,1,-1,1,-1,0)) 
FDR <- p.adjust(lrt24vs48 $table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrt24vs48 ))  
topTags(lrt24vs48 )
results <- topTags(lrt24vs48, n=1152)
write.table(as.matrix(results$table),file="DE_24vs48.txt",sep="\t")
dt.24vs48tags <- rownames(DGE[as.logical(dt)])
plotSmear(lrt24vs48 , de.tags=dt.24vs48tags)
abline(h=c(-2,2), col = "purple")

##Some time controlled comparisons (time 24)
lrtC24vsIi24 <- glmLRT(fitA, contrast = c(0,0,0,0,1,-1,0,0,0,0)) 
FDR <- p.adjust(lrtC24vsIi24 $table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrtC24vsIi24))  
topTags(lrtC24vsIi24)
results <- topTags(lrtC24vsIi24, n=69)
write.table(as.matrix(results$table),file="DE_C24vsIi24.txt",sep="\t")
dt.C24vsIi24tags <- rownames(DGE[as.logical(dt)])
plotSmear(lrtC24vsIi24, de.tags=dt.C24vsIi24tags)
abline(h=c(-2,2), col = "purple")

lrtC24vsIn24 <- glmLRT(fitA, contrast = c(0,0,0,0,1,0,0,-1,0,0)) 
FDR <- p.adjust(lrtC24vsIn24$table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrtC24vsIn24))  
results <- topTags(lrtC24vsIn24, n=48)
head(results)
write.table(as.matrix(results$table),file="DE_C24vsIn24.txt",sep="\t")
dt.C24vsIn24tags <- rownames(DGE[as.logical(dt)])
plotSmear(lrtC24vsIn24, de.tags=dt.C24vsIn24tags)
abline(h=c(-2,2), col = "purple")

lrtIi24vsIn24 <- glmLRT(fitA, contrast = c(0,0,0,0,0,1,0,-1,0,0)) 
FDR <- p.adjust(lrtIi24vsIn24$table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrtIi24vsIn24))  
topTags(lrtIi24vsIn24)
results <- topTags(lrtIi24vsIn24, n=1)
write.table(as.matrix(results$table),file="DE_Ii24vsIn24.txt",sep="\t")
dt.Ii24vsIn24tags <- rownames(DGE[as.logical(dt)])
plotSmear(lrtIi24vsIn24, de.tags=dt.Ii24vsIn24tags)
abline(h=c(-2,2), col = "purple")

lrtIi48vsIn48 <- glmLRT(fitA, contrast = c(0,0,0,0,0,0,1,0,-1,0)) 
FDR <- p.adjust(lrtIi48vsIn48$table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrtIi48vsIn48))  
topTags(lrtIi48vsIn48)
results <- topTags(lrtIi24vsIn24, n=1)

lrtW0vsC0 <- glmLRT(fitA, coef = 4) 
lrt2W0vsC0 <- glmLRT(fitA, contrast = c(0,0,0,-1,0,0,0,0,0,0)) 
FDR <- p.adjust(lrtW0vsC0$table$PValue, method = "BH")
FDR <- p.adjust(lrt2W0vsC0$table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrtW0vsC0))  
summary(dt <- decideTestsDGE(lrt2W0vsC0))  
topTags(lrtW0vsC0)
results <- topTags(lrtW0vsC0, n=23)
write.table(as.matrix(results$table),file="DE_W0vsC0.txt",sep="\t")
dt.W0vsC0tags <- rownames(DGE[as.logical(dt)])
plotSmear(lrtW0vsC0, de.tags=dt.W0vsC0tags)
abline(h=c(-2,2), col = "purple")


lrtI0Ii48 <- glmLRT(fitA, contrast = c(0,0,0,0,0,0,1,0,0,-1)) 
FDR <- p.adjust(lrtI0Ii48$table$PValue, method = "BH")
summary(dt <- decideTestsDGE(lrtI0Ii48))  
topTags(lrtI0Ii48)
results <- topTags(lrtI0Ii48, n=3267)
write.table(as.matrix(results$table),file="DE_I0Ii48.txt",sep="\t")
dt.I0Ii48tags <- rownames(DGE[as.logical(dt)])
plotSmear(lrtI0Ii48, de.tags=dt.I0Ii48tags)
abline(h=c(-2,2), col = "purple")

#HEATMAP with Euclidean Distance 
library(lattice)
rgb.palette <- colorRampPalette(c("green", "darkred"), space = "Lab")  #switch the colors if you like

#use the normalized counts from your DGEobject...= $counts
mycounts<-DGE$counts
designA <- model.matrix(~0+Block+Treat)

#Calculate Dispersions
A <- estimateGLMCommonDisp(A, designA, verbose = T)  #average over all genes, output: Disp = 0.06087, BCV = 0.2467 (went down compared to when I had no model)
A <- estimateGLMTrendedDisp(A, designA) # allows for a possible abudance trend
A <- estimateGLMTagwiseDisp(A, designA) #per tag (per gene?)

#Designs 
fitHM <- glmFit(DGE, designA)
mycounts.HM<-fitHM$counts
mydists.HM<-dist(t(mycounts.HM))  

# convert to matrix
mydistmatrix<-as.matrix(mydists.HM) ; # spearman correlation

#heatmap
levelplot(mydistmatrix, xlab="", ylab="",col.regions=rgb.palette(120), scales=list(x=list(rot=90)), main="Euclidean distance")


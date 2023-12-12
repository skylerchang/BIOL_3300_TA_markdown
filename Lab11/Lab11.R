# get the table of read counts by indicating the path to the file based on your system !!!!!!

readcounts <- read.table("~/Desktop/Lab10/RNA_seq_featurecounts.txt",header=TRUE)

# One of the requirements of the assay () slots is that the row.names 
# correspond to the gene IDs and the col.names to the sample names
row.names(readcounts)<-readcounts$Geneid

# In addition , we need to exclude all columns that do not contain read counts
readcounts<-readcounts[,-c(1:6)]

# give meaningful sample names - this can be achieved via numerous approaches
# the one shown here is the least generic and most error - prone one!
names(readcounts)<-c("Control_1","Control_2","Control_3","Mov10_1","Mov10_2","Mov10_3")


# make a data.frame with meta-data where row.names should match the individual
# sample names
sample_info<-data.frame(condition = gsub("_[0-9]+"," ",names(readcounts)),
                        row.names = names(readcounts))

# IF NEEDED, install DESeq2,which is not available via install.packages(),
# but through bioconductor
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)

# generate the DESeqDataSet
DESeq.ds<-DESeqDataSetFromMatrix(countData = readcounts,
                                 colData = sample_info,
                                 design = ~ condition)

# investigate different sample sizes
colSums(counts(DESeq.ds)) 

# calculate the size factor and add it to the data set
DESeq.ds<-estimateSizeFactors(DESeq.ds)
sizeFactors(DESeq.ds)

##################### Lab11 ########################
DESeq.ds <- estimateDispersions(DESeq.ds)
DESeq.ds <- nbinomWaldTest(DESeq.ds)
# Running the DE analysis
results(DESeq.ds)

# tells you which types of values can be extracted with results()
resultsNames(DESeq.ds) 

DGE.results <- results(DESeq.ds,
                       independentFiltering = TRUE,
                       alpha = 0.05)
# first line indicates which comparison was done for the log2FC
summary(DGE.results)

# the DESeqResult object can basically be handled like a data.frame
table(DGE.results$padj < 0.05)

#  A p-value of histogram 
hist(DGE.results$pvalue,
     col = "grey" , border = "white" , xlab = " " , ylab = " " ,
     main = "frequencies of p - values")

#  MA-plot
plotMA(DGE.results, alpha = 0.05, colNonSig = "gray60", colSig = "red",
       main = "Test: p.value < 0.05", ylim = c(-4,4))


# Heatmap plot
# load the library with the aheatmap () function
install.packages("pheatmap")
install.packages("NMF")
library(pheatmap)
library(NMF)
# aheatmap needs a matrix of values , e.g. , a matrix of DE genes with the
# transformed read counts for each replicate
# sort the results according to the adjusted p- value
DGE.results.sorted <- DGE.results[order(DGE.results$padj), ]

# identify genes with the desired adjusted p- value cut -off
DGEgenes <- rownames(subset(DGE.results.sorted , padj < 0.05))

# extract the normalized read counts for DE genes into a matrix
hm.mat_DGEgenes <- log.norm.counts[DGEgenes, ]

# plot the normalized read counts of DE genes sorted by the adjusted p- value
aheatmap(hm.mat_DGEgenes , Rowv = NA , Colv = NA)


# combine the heatmap with hierarchical clustering
pheatmap(hm.mat_DGEgenes,
           Rowv = TRUE , Colv = TRUE,                # add dendrograms to rows and column
           distfun = "euclidean" , hclustfun = "average")
# scale the read counts per gene to emphasize the sample-type-specific
# differences
pheatmap (hm.mat_DGEgenes ,
          Rowv = TRUE , Colv = TRUE ,
          distfun = "euclidean" , hclustfun = "average",
          scale = "row" ) 




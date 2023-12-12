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

# ALWAYS CHECK YOUR DATA AFTER YOU MANIPULATED IT!
str(readcounts)
head(readcounts,n=3) 


# make a data.frame with meta-data where row.names should match the individual
# sample names
sample_info<-data.frame(condition = gsub("_[0-9]+"," ",names(readcounts)),
                          row.names = names(readcounts))
sample_info


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

# test what counts () returns
counts(DESeq.ds)

# remove genes without any counts
DESeq.ds<-DESeq.ds[rowSums(counts(DESeq.ds))>0,]

# investigate different sample sizes
colSums(counts(DESeq.ds)) 


# calculate the size factor and add it to the data set
DESeq.ds<-estimateSizeFactors(DESeq.ds)
sizeFactors(DESeq.ds)

# if you check colData () again , you see that this now contains the sizeFactors
colData(DESeq.ds)


# counts () allows you to immediately retrieve the _ normalized _ read counts
counts.sf_normalized <- counts(DESeq.ds, normalized=TRUE)


# transform size - factor normalized read counts to log2 scale using a
# pseudocount of 1
log.norm.counts <- log2(counts.sf_normalized + 1)


# to plot the following two images next to each other
par(mfrow=c(1,1))
# first , boxplots of non - transformed read counts (one per sample )
boxplot(counts.sf_normalized, 
          main = "untransformed read counts", 
          ylab = "read counts")

# box plots of log2-transformed read counts
boxplot(log.norm.counts,
          main="log2-transformed read counts",
          ylab = "log2(read counts)" )

# obtain regularized log - transformed values
DESeq.rlog <- rlog(DESeq.ds,blind = TRUE)

# DESeq2 also offers a convenience function based on ggplot2 to do PCA directly on a DESeqDataSet:
#install.packages("ggplot2")
library(ggplot2)

# PCA 
P <- plotPCA(DESeq.rlog)
P <- P + theme_bw() + ggtitle("Rlog transformed counts")
print (P)


# cor () calculates the correlation between columns of a matrix
distance.m_rlog <- as.dist(1 - cor(log.norm.counts,method = "pearson"))

# plot () can directly interpret the output of hclust ()
plot(hclust(distance.m_rlog),
       labels = colnames(log.norm.counts),
       main = "rlog transformed read counts distance:Pearson correlation")



# load the library with the aheatmap () function
# install.packages("NMF")
library(NMF)

# DESeq2 uses the levels of the condition to determine the order of the
#comparison
str(colData(DESeq.ds)$condition)
# set control as the first-level-factor
colData(DESeq.ds)$condition <- relevel(colData(DESeq.ds)$condition, "Control ")
# Now, running the DGE analysis is very simple:
DESeq.ds <- DESeq(DESeq.ds)

# The results() function lets you extract the base means across samples, moderated log2 fold changes,
# standard errors, test statistics etc. for every gene.
DGE.results <- results(DESeq.ds, independentFiltering = TRUE , alpha = 0.05)
summary(DGE.results)

# the DESeqResult object can basically be handled like a data . frame
head(DGE.results)
table(DGE.results$padj < 0.05)
rownames(subset(DGE.results , padj < 0.05))

#  Histogram of p-value are present in a data set
hist(DGE.results$pvalue,
     col = "grey" , border = "white" , xlab = " " , ylab = " " ,
     main = "frequencies of p - values")

# MA plot
plotMA(DGE.results, alpha = 0.05 , main = "Control vs. Mov10 Over expressed" ,
       ylim = c(-4 , 4))

# Heatmap

# aheatmap needs a matrix of values , e.g. , a matrix of DE genes with the
# transformed read counts for each replicate
# sort the results according to the adjusted p- value
DGE.results.sorted <- DGE.results[order(DGE.results$padj), ]

# identify genes with the desired adjusted p- value cut -off
DGEgenes <- rownames(subset(DGE.results.sorted , padj < 0.05))

# extract the normalized read counts for DE genes into a matrix
hm.mat_DGEgenes <- log.norm.counts[DGEgenes, ]

# plot the normalized read counts of DE genes sorted by the adjusted p- value
aheatmap(hm.mat_DGEgenes, Rowv = NA , Colv = NA)

# combine the heatmap with hierarchical clustering
jpeg(file = "heatmap1.jpeg", width = 1024, height = 768)
aheatmap(hm.mat_DGEgenes, Rowv = TRUE, Colv = TRUE,               
         distfun = "euclidean" , hclustfun="average")

dev.off()

jpeg(file = "heatmap2.jpeg", width = 1024, height = 768)
aheatmap(hm.mat_DGEgenes ,
          Rowv = TRUE , Colv = TRUE ,
          distfun = "euclidean" , hclustfun = "average",
          scale = "row" ) 
dev.off()

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# list columns that can be retrieved from the annotation data base
columns(org.Hs.eg.db)

# make a batch retrieval for all DE genes
anno <- select(org.Hs.eg.db,
                 keys = DGEgenes,keytype ="SYMBOL",
                columns = c("GENENAME","GO"))


# check whether Mov10 pops up among the top upregulated genes
DGE.results.sorted_logFC<-DGE.results[order(DGE.results$log2FoldChange), ]
DGEgenes_logFC <- rownames (subset(DGE.results.sorted_logFC, padj < 0.05))
anno[match(DGEgenes_logFC, anno$SYMBOL), ]


BiocManager::install("clusterProfiler")
library(clusterProfiler)
## clusterProfiler requires a sorted vector where the values correspond
## to the measure used for the sorting 

DGE.results <- DGE.results[order(-1*DGE.results$log2FoldChange),]
genes_for_cp <- DGE.results$log2FoldChange
names(genes_for_cp) <- row.names(DGE.results)
## run the gene set enrichment analysis
gsea_kegg <- clusterProfiler::gseKEGG(geneList = genes_for_cp, organism = 'hsa',
                                      nPerm = 1000, minGSSize = 10,
                                      pvalueCutoff = 1, verbose = FALSE)







DGE.results.sorted_logFC<-DGE.results.sorted_logFC[order(-1*DGE.results.sorted_logFC$log2FoldChange),]
genes_for_cp<-DGE.results.sorted_logFC$log2FoldChange
eg<-bitr(row.names(DGE.results.sorted_logFC),fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(genes_for_cp) <- eg$ENTREZID
genes_for_cp<-na.omit(genes_for_cp)
genes_for_cp<-sort(genes_for_cp, decreasing = TRUE)
## run the gene set enrichment analysis
gse <- gseGO(geneList = genes_for_cp,OrgDb = organism,
                                      nPerm = 1000, minGSSize = 10,
                                      pvalueCutoff = 1, verbose = FALSE)




names(genes_for_cp) <- row.names(DGE.results.sorted_logFC)
OrgDb = "org.Hs.eg.db"
organism = "org.Hs.eg.db"
gse <- gseGO(geneList=genes_for_cp, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db",  
             pAdjustMethod = "none" )



cnetplot(gsea_kegg,
         showCategory = 2, colorEdge = TRUE, foldChange = genes_for_cp) +
         scale_colour_gradient2(name = "log2FC",
                         low = "navyblue", high = "red", mid = "white")

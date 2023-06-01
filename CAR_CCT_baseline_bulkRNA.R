
library(DESeq2)
library(pheatmap)
library(apeglm)
library(ggplot2)
library(stringr)
library(ggrepel)
# read in the metadata sample sheet
sampletable <- read.table("CAR_SampleSheet_DataSet1.txt", header=T, sep="\t")
rownames(sampletable) <- sampletable$SampleName
nrow(sampletable)
ncol(sampletable)
head(sampletable)
sampletable$Group <- as.factor(sampletable$Group)

#read in gene annotation file
tx2gene <- read.table("transcripts_to_genes.txt", 
                      sep="\t",
                      header=F)

#remove gene version number of ENSMBL ID, so that it is matching to the original count table
tx2gene[,2] <- sub("\\..*", "", tx2gene[,2])
head(tx2gene)

#read in gene counts file
count_matrix <- read.table("aligned2_RawCounts_DataSet1.txt", 
                           header=T, 
                           sep="\t", 
                           row.names=1)
se_star_matrix <- DESeqDataSetFromMatrix(countData = count_matrix,
                                         colData = sampletable,
                                         design = ~ Group)

nrow(se_star_matrix)
#Filter out low expression genes(optional)
se_star_matrix_filtered <- se_star_matrix[rowSums(counts(se_star_matrix)) > 10, ]
nrow(se_star_matrix_filtered)

#Fit statistical model
se_star_matrix_filtered <- DESeq(se_star_matrix_filtered)

#Normalize counts for future use
norm_counts <- counts(se_star_matrix_filtered, normalized = TRUE)

#Add gene symbols
norm_counts_symbols <- merge(unique(tx2gene[,2:3]),
                             data.frame(ID=rownames(norm_counts), norm_counts),
                             by=1, all=F)
write.table(norm_counts_symbols,
            "normalized_counts_DataSet1_filtered.txt", 
            quote=F, col.names=T, 
            row.names=F, sep="\t")

#vst transformation
vsd <- vst(se_star_matrix_filtered)

#Calculate between-sample distance matrix

sampleDistMatrix <- as.matrix(dist(t(assay(vsd))))
pheatmap(sampleDistMatrix)
plotPCA(object = vsd,
        intgroup = "Group")

# check results names: depends on what was modeled. Here is the "Group"
resultsNames(se_star_matrix_filtered)

#Extract results for AZ39d_R2 and AZ39_R2
de <- results(object =se_star_matrix_filtered, 
              name="Group_AZ39d_R2_vs_AZ39_R2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm")


de_shrink <- lfcShrink(dds = se_star_matrix_filtered,
                       coef="Group_AZ39d_R2_vs_AZ39_R2",
                       type="apeglm")
head(de)
head(de_shrink)

de_symbols <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(de_shrink), de_shrink), by=1, all=F)
write.table(de_symbols, "deseq2_DataSet1.txt", quote=F, col.names=T, row.names=F, sep="\t")

#Volcano

setwd("/Users/xiaoyuzhou/Desktop/aligned")

DataSet1 <- read.table("deseq2_DataSet1.txt", header = T)
colnames(DataSet1)[1:2] <- c("EnsmbleID", "Symbol")


DataSet1_filtered <- DataSet1[!str_detect(DataSet1$Symbol, "^TR[ABGD][VDJC]"), ]
nrow(DataSet1) - nrow(DataSet1_filtered) 
Removed_DataSet1_da <- DataSet1[str_detect(DataSet1$Symbol, "^TR[ABGD][VDJC]"), ]

DataSet1_filtered$Category <- "Not Sig"
DataSet1_filtered$Category[DataSet1_filtered$log2FoldChange > 1 & DataSet1_filtered$padj < 0.01] <- "Up"
DataSet1_filtered$Category[DataSet1_filtered$log2FoldChange < -1 & DataSet1_filtered$padj < 0.01] <- "Down"
DataSet1_da_up <- DataSet1_filtered[DataSet1_filtered $Category == "Up", ]
DataSet1_da_down <- DataSet1_filtered[DataSet1_filtered $Category == "Down", ]
DataSet1_da_de <- rbind(DataSet1_da_up, DataSet1_da_down)
head(DataSet1_da_de)
write.table(DataSet1_da_de,
            "DataSet1_DE_all.txt", 
            quote=F, col.names=T, 
            row.names=F, sep="\t")

DataSet1_filtered$delabel <- NA
DataSet1_filtered$delabel[DataSet1_filtered$Category != "Not Sig"] <- 
  DataSet1_filtered$Symbol[DataSet1_filtered$Category != "Not Sig"]


DataSet1_plot_da <- 
  ggplot(data=DataSet1_filtered, aes(x=log2FoldChange,y=-log10(padj), col=Category, label = delabel))+
  geom_point()+
  geom_text_repel()+
  scale_color_manual(values=c("cornflowerblue", "green3", "salmon"))


####heatmap####

#load normalized count table
DataSet1_norm_filtered <- read.table("normalized_counts_DataSet1_filtered.txt", header = T)
head(DataSet1_norm_filtered)
#annotate colnames
colnames(DataSet1_norm_filtered)[1:8] <- c("EnsmbleID", "Symbol", "39A2_1", "39A2_2", "39A2_3",
                                            "39D2_1", "39D2_2", "39D2_3")
#select only differentially expressed genes
all_de_uniq_NormCounts <- DataSet1_norm_filtered[DataSet1_norm_filtered$Symbol %in% DataSet1_da_de$Symbol, ]
all_de_uniq_NormCounts <- all_de_uniq_NormCounts[, -1]

#make a numeric matrix
rownames(all_de_uniq_NormCounts) <- all_de_uniq_NormCounts$Symbol
all_de_uniq_NormCounts <- all_de_uniq_NormCounts[, -1]
head(all_de_uniq_NormCounts)
all_de_uniq_NormCounts_matrix <- as.matrix(all_de_uniq_NormCounts)

#plot with pheatmap
DataSet1_heatmap <- pheatmap::pheatmap(all_de_uniq_NormCounts_matrix, scale="row", cluster_cols = F, 
                         clustering_method = 'ward.D2', fontsize_row = 0.5,cutree_rows = 2)
?pheatmap

#### DAVID ####

DAVID_UP <- read.table("/Users/xiaoyuzhou/Desktop/aligned/DAVID analysis/DataSet1_UP.txt",header = T, sep = "\t")
DAVID_Down <- read.table("/Users/xiaoyuzhou/Desktop/aligned/DAVID analysis/DataSet1_Down.txt",header = T, sep = "\t")
DAVID_UP_TOP <- DAVID_UP[1:6, ]
DAVID_Down_TOP <- DAVID_Down[1:6, ]


ggplot(data=DAVID_UP_TOP, aes(x=-log10(PValue),y= reorder(Term, -log10(PValue))))+
  geom_bar(stat = "identity", color = "black", fill="darkred", width = 0.5)+
  theme_minimal()

ggplot(data=DAVID_Down_TOP, aes(x=-log10(PValue),y= reorder(Term, -log10(PValue))))+
  geom_bar(stat = "identity", color = "black", fill="steelblue", width = 0.5)+
  theme_minimal()


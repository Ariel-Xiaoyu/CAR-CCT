
library(DESeq2)
library(pheatmap)
library(apeglm)
library(ggplot2)
library(stringr)
library(gplots)
library(RColorBrewer)
library(NMF)
library(ggrepel)
library(enrichplot)

# read in the metadata sample sheet
sampletable <- read.table("CAR_SampleSheet_DataSet2.txt", header=T, sep="\t")
rownames(sampletable) <- sampletable$SampleName
nrow(sampletable)
ncol(sampletable)
head(sampletable)
sampletable$Group <- as.factor(sampletable$Group)

#read in gene annotation file
tx2gene <- read.table("transcripts_to_genes.txt", 
                      sep="\t",
                      header=F)
tx2gene[,2] <- sub("\\..*", "", tx2gene[,2])
head(tx2gene)
#read in gene counts file
count_matrix <- read.table("aligned2_RawCounts_DataSet2.txt", 
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
#norm_counts <- log2(counts(se_star_matrix, normalized = TRUE)+1)
norm_counts_filtered_wt <- counts(se_star_matrix_filtered, normalized = TRUE)

#Add gene symbols
norm_counts_symbols <- merge(unique(tx2gene[,2:3]),
                             data.frame(ID=rownames(norm_counts_filtered_wt), norm_counts_filtered_wt),
                             by=1, all=F)
write.table(norm_counts_symbols,
            "normalized_counts_DataSet2_filtered_wt.txt", 
            quote=F, col.names=T, 
            row.names=F, sep="\t")

#vst transformation
vsd <- vst(se_star_matrix_filtered)

#Calculate between-sample distance matrix

sampleDistMatrix <- as.matrix(dist(t(assay(vsd))))
pheatmap(sampleDistMatrix)
PCA_basline <- plotPCA(object = vsd,
        intgroup = "Group")
ggsave(plot = PCA_basline, filename = "PCA_basline.pdf", 
       dpi = 300, width = 8, height = 4, units = 'in')

# check results names: depends on what was modeled. Here is the "Group"
resultsNames(se_star_matrix_filtered)

#Extract results for AZ39d_R2 and AZ39_R2
#de <- results(object =se_star_matrix, name="Group_AZ39d_R2_vs_AZ39_R2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm")


de_shrink_ba <- lfcShrink(dds = se_star_matrix_filtered,
                       coef="Group_AZ39b_vs_AZ39",
                       type="apeglm")
de_shrink_ca <- lfcShrink(dds = se_star_matrix_filtered,
                          coef="Group_AZ39c_vs_AZ39",
                          type="apeglm")
de_shrink_da <- lfcShrink(dds = se_star_matrix_filtered,
                          coef="Group_AZ39d_vs_AZ39",
                          type="apeglm")


de_symbols_ba <- merge(unique(tx2gene[,2:3]),
                       data.frame(ID=rownames(de_shrink_ba), de_shrink_ba), 
                       by=1, all=F)
de_symbols_ca <- merge(unique(tx2gene[,2:3]),
                       data.frame(ID=rownames(de_shrink_ca), de_shrink_ca), 
                       by=1, all=F)
de_symbols_da <- merge(unique(tx2gene[,2:3]),
                       data.frame(ID=rownames(de_shrink_da), de_shrink_da), 
                       by=1, all=F)


write.table(de_symbols_ba, "deseq2_DataSet2_AZ39b_vs_AZ39.txt", 
            quote=F, col.names=T, row.names=F, sep="\t")
write.table(de_symbols_ca, "deseq2_DataSet2_AZ39c_vs_AZ39.txt", 
            quote=F, col.names=T, row.names=F, sep="\t")
write.table(de_symbols_da, "deseq2_DataSet2_AZ39d_vs_AZ39.txt", 
            quote=F, col.names=T, row.names=F, sep="\t")

#Volcano
DataSet2_ba <- read.table("deseq2_DataSet2_AZ39b_vs_AZ39.txt", header = T)
DataSet2_ca <- read.table("deseq2_DataSet2_AZ39c_vs_AZ39.txt", header = T)
DataSet2_da <- read.table("deseq2_DataSet2_AZ39d_vs_AZ39.txt", header = T)

colnames(DataSet2_ba)[1:2] <- c("EnsmbleID", "Symbol")
colnames(DataSet2_ca)[1:2] <- c("EnsmbleID", "Symbol")
colnames(DataSet2_da)[1:2] <- c("EnsmbleID", "Symbol")

# Set DE cutoff
DataSet2_ba $Category <- "Not Sig"
DataSet2_ba $Category[DataSet2_ba$log2FoldChange > 0.5 & DataSet2_ba$padj < 0.01] <- "Up"
DataSet2_ba$Category[DataSet2_ba$log2FoldChange < -0.5 & DataSet2_ba$padj < 0.01] <- "Down"
ba_up <- DataSet2_ba[DataSet2_ba $Category == "Up", ]
ba_down <- DataSet2_ba[DataSet2_ba $Category == "Down", ]
ba_de <- rbind(ba_up, ba_down)

DataSet2_ca $Category <- "Not Sig"
DataSet2_ca $Category[DataSet2_ca$log2FoldChange > 0.5 & DataSet2_ca$padj < 0.01] <- "Up"
DataSet2_ca$Category[DataSet2_ca$log2FoldChange < -0.5 & DataSet2_ca$padj < 0.01] <- "Down"
ca_up <- DataSet2_ca[DataSet2_ca $Category == "Up", ]
ca_down <- DataSet2_ca[DataSet2_ca $Category == "Down", ]
ca_de <- rbind(ca_up, ca_down)

DataSet2_da $Category <- "Not Sig"
DataSet2_da $Category[DataSet2_da$log2FoldChange > 0.5 & DataSet2_da$padj < 0.01] <- "Up"
DataSet2_da$Category[DataSet2_da$log2FoldChange < -0.5 & DataSet2_da$padj < 0.01] <- "Down"
da_up <- DataSet2_da[DataSet2_da $Category == "Up", ]
da_down <- DataSet2_da[DataSet2_da $Category == "Down", ]
da_de <- rbind(da_up, da_down)

DataSet2_ba$delabel <- NA
DataSet2_ba$delabel[DataSet2_ba$Category != "Not Sig"] <- 
  DataSet2_ba$Symbol[DataSet2_ba$Category != "Not Sig"]

DataSet2_ca$delabel <- NA
DataSet2_ca$delabel[DataSet2_ca$Category != "Not Sig"] <- 
  DataSet2_ca$Symbol[DataSet2_ca$Category != "Not Sig"]

DataSet2_da$delabel <- NA
DataSet2_da$delabel[DataSet2_da$Category != "Not Sig"] <- 
  DataSet2_da$Symbol[DataSet2_da$Category != "Not Sig"]
all_de_df <- rbind(ba_de, ca_de, da_de)

CAR1_overlap_bc_down <- intersect(ba_down$Symbol, ca_down$Symbol)
CAR1_overlab_bcd_down <- intersect(CAR1_overlap_bc, da_down$Symbol)


write.table(all_de_df,
            "DataSet2_DE_all_woTCRfilter.txt", 
            quote=F, col.names=T, 
            row.names=F, sep="\t")
####heatmap####

#select DE list
all_de <- c(ba_de$Symbol, ca_de$Symbol, da_de$Symbol)
all_de_uniq <- unique(all_de)

#load normalized count table
DataSet2_norm_filtered <- read.table("normalized_counts_DataSet2_filtered_wt.txt", header = T)

#annotate colnames
colnames(DataSet2_norm_filtered)[1:14] <- c("EnsmbleID", "Symbol", "39A0_1", "39A0_2", "39A0_3",
                                            "39B0_1", "39B0_2", "39B0_3", "39C0_1", "39C0_2", "39C0_3",
                                            "39D0_1", "39D0_2", "39D0_3")
#select only differentially expressed genes
all_de_uniq_NormCounts <- DataSet2_norm_filtered[DataSet2_norm_filtered$Symbol %in% all_de_uniq, ]
all_de_uniq_NormCounts <- all_de_uniq_NormCounts[, -1]

#make a numeric matrix
rownames(all_de_uniq_NormCounts) <- all_de_uniq_NormCounts$Symbol
all_de_uniq_NormCounts <- all_de_uniq_NormCounts[, -1]
head(all_de_uniq_NormCounts)

#plot with gplot


all_de_uniq_NormCounts_matrix <- as.matrix(all_de_uniq_NormCounts)
heatmap.2(all_de_uniq_NormCounts_matrix, col=rev(brewer.pal(9,"RdBu")), scale="row", Colv = "NA")

#plot with pheatmap
p1 <- pheatmap::pheatmap(all_de_uniq_NormCounts_matrix, scale="row", cluster_cols = F, 
                   clustering_method = 'ward.D2', fontsize_row = 4,cutree_rows = 4)
geneModule <- list()
for(i in 1:4) geneModule[[i]] <- names(cutree(p1$tree_row, k = 4))[which(cutree(p1$tree_row, k = 4) == i)]
for(i in 1:4) write.table(geneModule[[i]], file = paste0('/Users/xiaoyuzhou/Desktop/aligned/GeneModule_',i,'.txt'), 
                          sep = '\t', col.names = F,row.names = F, quote = F)

#plot with NMF

NMF_heatmap <- aheatmap(all_de_uniq_NormCounts_matrix, color = rev(brewer.pal(9,"RdBu")), 
         scale="row", annColors = "Set2", annRow=all_de_uniq_NormCounts$Symbol, 
         Colv = NA, filename = "NMF_heatmap.pdf")

#find shared downregulated genes in each pair-wised comparison 
bc_overlap <- intersect(ba_down$Symbol, ca_down$Symbol)
cd_overlap <- intersect(ca_down$Symbol, da_down$Symbol)
bd_overlap <- intersect(ba_down$Symbol, da_down$Symbol)
bcd_overlap <- intersect(bc_overlap, da_down$Symbol)

####Volcano####

ggplot(data=DataSet2_ba, aes(x=log2FoldChange,y=-log10(padj), col=Category, label = delabel))+
  geom_point()+
  geom_text_repel()+
  scale_color_manual(values=c("cornflowerblue", "green3", "salmon"))


ggplot(data=DataSet2_ca, aes(x=log2FoldChange,y=-log10(padj), col=Category, label = delabel))+
  geom_point()+
  geom_text_repel()+
  scale_color_manual(values=c("cornflowerblue", "green3", "salmon"))

ggplot(data=DataSet2_da, aes(x=log2FoldChange,y=-log10(padj), col=Category, label = delabel))+
  geom_point()+
  geom_text_repel()+
  scale_color_manual(values=c("cornflowerblue", "green3", "salmon"))

#####pathway analysis####

#bcd_all_down
CAR1_all_down <- rbind(ba_down, ca_down, da_down)
CAR1_bc_down <- rbind(ba_down, ca_down)

CAR1_all_down_ENTREZID <- select(org.Hs.eg.db, keys = CAR1_all_down$Symbol,
                         columns = c("ENTREZID", "SYMBOL"),
                         keytype = "SYMBOL")

CAR1_all_down_go <- data.frame(CAR1_all_down_ENTREZID$ENTREZID, 
                               CAR1_all_down$log2FoldChange)
CAR1_all_down_go <- CAR1_all_down_go[!is.na(CAR1_all_down_go[,1]),]
geneList_CAR1_all_down <- CAR1_all_down_go[,2]
names(geneList_CAR1_all_down) <- as.character(CAR1_all_down_go[,1])
geneList_CAR1_all_down <- sort(geneList_CAR1_all_down, decreasing = TRUE)
class(geneList_CAR1_all_down)


ego_CAR1_all_down <- enrichGO(names(geneList_CAR1_all_down), 
                      OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
x <- data.frame(ego_CAR1_all_down)
write.table(x, "CAR1_down_all_gesea.txt", sep = "\t", quote = F)

#select immune related categories
ego_CAR1_all_down_select <- 
  ego_CAR1_all_down[c(155, 162, 164, 172, 177, 180, 184, 216),asis=T]

CAR1_all_down_gsea <- dotplot(ego_CAR1_all_down_select, showCategory= 8) + ggtitle("dotplot for ORA")
ggsave(plot = CAR1_all_down_gsea, filename = "CAR1_all_down_gsea.pdf", 
       dpi = 300, width = 8, height = 8, units = 'in')

#select top 10 catagories
dotplot(ego_CAR1_all_down, showCategory=10) + ggtitle("dotplot for ORA")


#bc_down
CAR1_bc_down_ENTREZID <-  select(org.Hs.eg.db, keys = CAR1_bc_down$Symbol,
                                 columns = c("ENTREZID", "SYMBOL"),
                                 keytype = "SYMBOL")

CAR1_bc_down_go <- data.frame(CAR1_bc_down_ENTREZID$ENTREZID, 
                              CAR1_bc_down$log2FoldChange)
CAR1_bc_down_go <- CAR1_bc_down_go[!is.na(CAR1_bc_down_go[,1]),]
geneList_CAR1_bc_down_go <- CAR1_bc_down_go[,2]
names(geneList_CAR1_bc_down_go) <- as.character(CAR1_bc_down_go[,1])
geneList_CAR1_bc_down_go <- sort(geneList_CAR1_bc_down_go, decreasing = TRUE)
class(geneList_CAR1_bc_down_go)

ego_CAR1_bc_down <- enrichGO(names(geneList_CAR1_bc_down_go), 
                              OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
y <- data.frame(ego_CAR1_bc_down)
write.table(y, "CAR1_bc_down_gesea.txt", sep = "\t", quote = F)



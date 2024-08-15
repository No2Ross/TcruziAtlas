source("own_method_functions.R")
#Trade seq and STACAS do not get along 
library(phateR)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(phateR)
library(dplyr)
library(Seurat)
library(rgl)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(viridis)
library(scran)
library(RColorBrewer)
library(SingleCellExperiment)
library(reshape2)
library(slingshot)
library(stats)
library(stringi)
library(harmony)
library(stringr)
library(mclust)
library(gplots)
library(cvms)
library(candisc)
library(EnhancedVolcano)

library(unixtools)

setwd("/datastore/Ross")

#set.tempdir("pooling-scRNA/")
ulimit::memory_limit(50000)

setwd("/datastore/Ross")

de_dir <- "/datastore/Ross/cruzi_files/lifecycle/bulk_RNA_seq/output/DE_tables/"


names_objects <- c("EarlyAmast1", "EarlyAmast2", "MidAmast1", "MidAmast2", "LateAmast1", "LateAmast2",
                   "Trypomast1", "Trypomast3", "Trypomast4", "Epimast1", "Epimast2","Metacyclic1","Metacyclic2")

first <- read.delim("/datastore/Ross/cruzi_paper/bulk/input_objects/marta_data/Epimast1_counts.txt", sep ="")
first <- first[,1:7]
colnames(first) <- first[1,]
first <- first[-1,]
first$Geneid <- str_replace_all(first$Geneid, "TcruziDm28cUTR_MT_", "")
row.names(first) <- first$Geneid
first <- first[,1]



exp_mtx <- matrix(0, ncol = length(names_objects), nrow = length(first))

for(i in 1:length(names_objects)){
  x <- read.delim(paste0("/datastore/Ross/cruzi_paper/bulk/input_objects/marta_data/", names_objects[i], "_counts.txt"), sep ="")
  x <- x[,1:7]
  colnames(x) <- x[1,]
  x <- x[-1,]
  print(names_objects[i])
  x$Geneid <- str_replace_all(x$Geneid, "TcruziDm28cUTR_MT_", "")
  row.names(x) <- x$Geneid
  
  exp_mtx[,i] <- as.numeric(x[,7])
}

colnames(exp_mtx) <- names_objects
row.names(exp_mtx) <- first

exp_mtx <- exp_mtx[ str_detect(row.names(exp_mtx), "Hsampien", negate = TRUE), ]

exp_mtx <- exp_mtx[-c(15320,15321),]

metaData_obj <- data.frame("stage_fine" = str_replace_all(names_objects, "[1234]", ""),
                           "stage_coarse" = as.factor( c("Amast", "Amast", "Amast", "Amast", "Amast", "Amast",
                                                         "Trypomast", "Trypomast", "Trypomast", "Epimast", "Epimast",
                                                         "Metacyclic", "Metacyclic")), 
                           "sample" =  as.factor( c("EarlyAmast1", "EarlyAmast2", "MidAmast1", "MidAmast2", "LateAmast1", "LateAmast2",
                                                    "Trypomast1", "Trypomast2", "Trypomast2", 
                                                    "Epimast1", "Epimast2","Metacyclic1","Metacyclic2")),
                           
                           "sample_replicate" = as.factor( c("EarlyAmast1", "EarlyAmast2", "MidAmast1", "MidAmast2", "LateAmast1", "LateAmast2",
                                                             "Trypomast1",  "Trypomast2_1", "Trypomast2_2", 
                                                             "Epimast1", "Epimast2","Metacyclic1","Metacyclic2")),
                           
                           "shortSample" = as.factor( c("eA1", "eA2", "mA1", "mA2", "lA1", "lA2",
                                                        "T1","T2","T2","E1", "E2", "M1", "M2")))



exp_mtx <- exp_mtx[,c("EarlyAmast1", "EarlyAmast2", "MidAmast1", "MidAmast2", "LateAmast1", "LateAmast2","Trypomast1",
                      "Trypomast3","Trypomast4","Epimast1", "Epimast2","Metacyclic1","Metacyclic2")]

metaData_obj$stage_fine <- factor(metaData_obj$stage_fine, levels = c("EarlyAmast", "MidAmast", "LateAmast",
                                                                      "Trypomast", "Epimast", "Metacyclic"))

metaData_obj$stage_coarse <- factor(metaData_obj$stage_coarse, levels = c("Amast", "Trypomast",
                                                                          "Epimast", "Metacyclic"))

metaData_obj$sample <- factor(metaData_obj$sample, levels = c("EarlyAmast1", "EarlyAmast2", "MidAmast1", "MidAmast2", "LateAmast1", "LateAmast2",
                                                              "Trypomast1","Trypomast2", "Epimast1", "Epimast2", "Metacyclic1", "Metacyclic2"))

metaData_obj$sample_replicate <- factor(metaData_obj$sample_replicate, levels = c("EarlyAmast1", "EarlyAmast2", "MidAmast1", "MidAmast2", "LateAmast1", "LateAmast2",
                                                                                  "Trypomast1", "Trypomast2_1","Trypomast2_2", "Epimast1", "Epimast2", "Metacyclic1", "Metacyclic2"))


metaData_obj$shortSample <- factor(metaData_obj$shortSample, levels = c("eA1", "eA2", "mA1", "mA2", "lA1", "lA2","T1", "T2","E1", "E2", "M1", "M2"))


library(DESeq2)  

dds_cruzi <- DESeqDataSetFromMatrix(exp_mtx,
                                    colData = metaData_obj,
                                    design = ~0 + stage_coarse)

dds_cruzi <- collapseReplicates(dds_cruzi, dds_cruzi$sample)

smallestGroupSize <- 11
keep <- rowSums(counts(dds_cruzi) >= 10) >= smallestGroupSize
dds_cruzi <- dds_cruzi[keep,]

dds_cruzi <- DESeq(dds_cruzi)
res <- results(dds_cruzi)

vsdata <- vst(dds_cruzi, blind=FALSE)

pdf("/datastore/Ross/cruzi_paper/bulk/paper/plots/BulkRNA_deseq2_PCA_stages.pdf")
plotPCA(vsdata, intgroup="stage_fine")
dev.off()

pdf("/datastore/Ross/cruzi_paper/bulk/paper/plots/BulkRNA_deseq2_PCA_sample.pdf")
plotPCA(vsdata, intgroup="sample")
dev.off()

plotPCA(vsdata, intgroup="stage_fine")

#stage Coarse contrasts

metacyclic_markers <- results(dds_cruzi, 
                              contrast=list(c("stage_coarseMetacyclic"), 
                                            c("stage_coarseEpimast","stage_coarseTrypomast","stage_coarseAmast")),
                              listValues=c(1, -1/3))

metacyclic_markers <- subset(metacyclic_markers, log2FoldChange > 0.5 & padj < 0.05)

amast_markers <- results(dds_cruzi, 
                              contrast=list(c("stage_coarseAmast"), 
                                            c("stage_coarseEpimast","stage_coarseTrypomast","stage_coarseMetacyclic")),
                              listValues=c(1, -1/3))

amast_markers <- subset(amast_markers, log2FoldChange > 0.5 & padj < 0.05)


epimast_markers <- results(dds_cruzi, 
                         contrast=list(c("stage_coarseEpimast"), 
                                       c("stage_coarseAmast","stage_coarseTrypomast","stage_coarseMetacyclic")),
                         listValues=c(1, -1/3))

epimast_markers <- subset(epimast_markers, log2FoldChange > 0.5 & padj < 0.05)


trypomast_markers <- results(dds_cruzi, 
                           contrast=list(c("stage_coarseTrypomast"), 
                                         c("stage_coarseAmast","stage_coarseEpimast","stage_coarseMetacyclic")),
                           listValues=c(1, -1/3))

trypomast_markers <- subset(trypomast_markers, log2FoldChange > 0.5 & padj < 0.05)

#Remove non-unique genes
trypomast_genes <- setdiff(row.names(trypomast_markers), unique( c( row.names(metacyclic_markers),
                                                            row.names(amast_markers),
                                                            row.names(epimast_markers) ) ) )

epimast_genes <- setdiff(row.names(epimast_markers), unique( c( row.names(metacyclic_markers),
                                                                    row.names(amast_markers),
                                                                    row.names(trypomast_markers) ) ) )

amast_genes <- setdiff(row.names(amast_markers), unique( c( row.names(metacyclic_markers),
                                                                row.names(epimast_markers),
                                                                row.names(trypomast_markers) ) ) )

meta_genes <- setdiff(row.names(metacyclic_markers), unique( c( row.names(epimast_markers),
                                                                row.names(amast_markers),
                                                                row.names(trypomast_markers) ) ) )

metacyclic_markers <- metacyclic_markers[meta_genes,]
amast_markers <- amast_markers[amast_genes,]
epimast_markers <- epimast_markers[epimast_genes,]
trypomast_markers <- trypomast_markers[trypomast_genes,]

write.csv(metacyclic_markers,
          "/datastore/Ross/cruzi_paper/bulk/paper/DEG_files/stageCoarse/metacyclics_unique_sigGenes.csv")

write.csv(trypomast_markers,
          "/datastore/Ross/cruzi_paper/bulk/paper/DEG_files/stageCoarse/trypomast_unique_sigGenes.csv")

write.csv(epimast_markers,
          "/datastore/Ross/cruzi_paper/bulk/paper/DEG_files/stageCoarse/epimast_unique_sigGenes.csv")

write.csv(amast_markers,
          "/datastore/Ross/cruzi_paper/bulk/paper/DEG_files/stageCoarse/amast_unique_sigGenes.csv")

dim(metacyclic_markers)
dim(trypomast_markers)
dim(epimast_markers)
dim(amast_markers)

plotCounts(dds_cruzi, gene="C4B63_135g21", intgroup="stage_coarse") 
plotCounts(dds_cruzi, gene="C4B63_9g448", intgroup="stage_coarse") 

amast_df <- as.data.frame(amast_markers)
amast_logfc_order <- amast_df[order(amast_df$log2FoldChange, decreasing = T),]

epimast_df <- as.data.frame(epimast_markers)
epimast_logfc_order <- epimast_df[order(epimast_df$log2FoldChange, decreasing = T),]

trypomast_df <- as.data.frame(trypomast_markers)
trypomast_logfc_order <- trypomast_df[order(trypomast_df$log2FoldChange, decreasing = T),]

meta_df <- as.data.frame(metacyclic_markers)
meta_logfc_order <- meta_df[order(meta_df$log2FoldChange, decreasing = T),]

#Heatmap
sigGenes <- unique(c(row.names(amast_logfc_order)[1:10], row.names(trypomast_logfc_order)[1:10],  
                     row.names(epimast_logfc_order)[1:10], row.names(meta_logfc_order)[1:10]))


ntd <- normTransform(dds_cruzi)
library("pheatmap")
df <- as.data.frame(colData(dds_cruzi)[,c("stage_fine","stage_coarse")])

mtx_test <- as.matrix(assay(ntd)[sigGenes,])

row.subsections <- c(10,10,10,10)
row_split = data.frame("order"= rep(c("Amastigote markers", "Trypomastigote markers","Epimastigote markers", "Metacyclic markers"), row.subsections))
row_split$order <- factor(row_split$order, c("Amastigote markers", "Trypomastigote markers","Epimastigote markers", "Metacyclic markers"))

ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows = F,
                        row_split = row_split, cluster_row_slices = F)

mtx_test <- t(scale(t(as.matrix(assay(ntd)[sigGenes,]))))

pdf("/datastore/Ross/cruzi_paper/bulk/paper/plots/BulkRNA_deseq2_top10_coarse_lifecycleMarkers_scaledHeatmap.pdf",
    height = 10)
ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows = F,
                        row_split = row_split, cluster_row_slices = F,
                        heatmap_legend_param = list(title = "Scaled expression"))
dev.off()



#Stage fine comparison

dds_cruzi <- DESeqDataSetFromMatrix(exp_mtx,
                                    colData = metaData_obj,
                                    design = ~0 + stage_fine)

dds_cruzi <- collapseReplicates(dds_cruzi, dds_cruzi$sample)

smallestGroupSize <- 11
keep <- rowSums(counts(dds_cruzi) >= 10) >= smallestGroupSize
dds_cruzi <- dds_cruzi[keep,]

dds_cruzi <- DESeq(dds_cruzi)
res <- results(dds_cruzi)

vsdata <- vst(dds_cruzi, blind=FALSE)

# pdf("/datastore/Ross/cruzi_paper/bulk/plots/BulkRNA_deseq2_PCA_stages.pdf")
# plotPCA(vsdata, intgroup="stage_fine")
# dev.off()
# 
# pdf("/datastore/Ross/cruzi_paper/bulk/plots/BulkRNA_deseq2_PCA_sample.pdf")
# plotPCA(vsdata, intgroup="sample")
# dev.off()

plotPCA(vsdata, intgroup="stage_fine")

metacyclic_markers <- results(dds_cruzi, 
                              contrast=list(c("stage_fineMetacyclic"), 
                                            c("stage_fineEpimast","stage_fineTrypomast",
                                              "stage_fineEarlyAmast", "stage_fineMidAmast", "stage_fineLateAmast")),
                              listValues=c(1, -1/5))

metacyclic_markers <- subset(metacyclic_markers, log2FoldChange > 0.5 & padj < 0.05)

amastE_markers <- results(dds_cruzi, 
                          contrast=list(c("stage_fineEarlyAmast"), 
                                        c("stage_fineEpimast","stage_fineTrypomast",
                                          "stage_fineMetacyclic", "stage_fineMidAmast", "stage_fineLateAmast")),
                          listValues=c(1, -1/5))

amastE_markers <- subset(amastE_markers, log2FoldChange > 0.5 & padj < 0.05)

amastM_markers <- results(dds_cruzi, 
                          contrast=list(c("stage_fineMidAmast"), 
                                        c("stage_fineEpimast","stage_fineTrypomast",
                                          "stage_fineEarlyAmast", "stage_fineMetacyclic", "stage_fineLateAmast")),
                          listValues=c(1, -1/5))

amastM_markers <- subset(amastM_markers, log2FoldChange > 0.5 & padj < 0.05)


amastL_markers <- results(dds_cruzi, 
                          contrast=list(c("stage_fineLateAmast"), 
                                        c("stage_fineEpimast","stage_fineTrypomast",
                                          "stage_fineEarlyAmast", "stage_fineMidAmast", "stage_fineMetacyclic")),
                          listValues=c(1, -1/5))

amastL_markers <- subset(amastL_markers, log2FoldChange > 0.5 & padj < 0.05)


epimast_markers <- results(dds_cruzi, 
                           contrast=list(c("stage_fineEpimast"), 
                                         c("stage_fineMetacyclic","stage_fineTrypomast",
                                           "stage_fineEarlyAmast", "stage_fineMidAmast", "stage_fineLateAmast")),
                           listValues=c(1, -1/5))

epimast_markers <- subset(epimast_markers, log2FoldChange > 0.5 & padj < 0.05)


trypomast_markers <- results(dds_cruzi, 
                             contrast=list(c("stage_fineTrypomast"), 
                                           c("stage_fineEpimast","stage_fineMetacyclic",
                                             "stage_fineEarlyAmast", "stage_fineMidAmast", "stage_fineLateAmast")),
                             listValues=c(1, -1/5))

trypomast_markers <- subset(trypomast_markers, log2FoldChange > 0.5 & padj < 0.05)

#Remove non-unique genes
trypomast_genes <- setdiff(row.names(trypomast_markers), unique( c( row.names(metacyclic_markers),
                                                                    row.names(amastE_markers),
                                                                    row.names(amastM_markers),
                                                                    row.names(amastL_markers),
                                                                    row.names(epimast_markers) ) ) )

epimast_genes <- setdiff(row.names(epimast_markers), unique( c( row.names(metacyclic_markers),
                                                                row.names(amastE_markers),
                                                                row.names(amastM_markers),
                                                                row.names(amastL_markers),
                                                                row.names(trypomast_markers) ) ) )

amastE_genes <- setdiff(row.names(amastE_markers), unique( c( row.names(metacyclic_markers),
                                                            row.names(trypomast_markers),
                                                            row.names(amastM_markers),
                                                            row.names(amastL_markers),
                                                            row.names(epimast_markers) ) ) )

amastM_genes <- setdiff(row.names(amastM_markers), unique( c( row.names(metacyclic_markers),
                                                              row.names(amastE_markers),
                                                              row.names(trypomast_markers),
                                                              row.names(amastL_markers),
                                                              row.names(epimast_markers) ) ) )

amastL_genes <- setdiff(row.names(amastL_markers), unique( c( row.names(metacyclic_markers),
                                                              row.names(amastE_markers),
                                                              row.names(amastM_markers),
                                                              row.names(trypomast_markers),
                                                              row.names(epimast_markers) ) ) )

meta_genes <- setdiff(row.names(metacyclic_markers), unique( c( row.names(trypomast_markers),
                                                                row.names(amastE_markers),
                                                                row.names(amastM_markers),
                                                                row.names(amastL_markers),
                                                                row.names(epimast_markers) ) ) )

metacyclic_markers <- metacyclic_markers[meta_genes,]
amastE_markers <- amastE_markers[amastE_genes,]
amastM_markers <- amastM_markers[amastM_genes,]
amastL_markers <- amastL_markers[amastL_genes,]
epimast_markers <- epimast_markers[epimast_genes,]
trypomast_markers <- trypomast_markers[trypomast_genes,]

write.csv(metacyclic_markers,
          "/datastore/Ross/cruzi_paper/bulk/paper/DEG_files/stageFine/metacyclics_unique_sigGenes.csv")

write.csv(trypomast_markers,
          "/datastore/Ross/cruzi_paper/bulk/paper/DEG_files/stageFine/trypomast_unique_sigGenes.csv")

write.csv(epimast_markers,
          "/datastore/Ross/cruzi_paper/bulk/paper/DEG_files/stageFine/epimast_unique_sigGenes.csv")

write.csv(amastE_markers,
          "/datastore/Ross/cruzi_paper/bulk/paper/DEG_files/stageFine/amastEarly_unique_sigGenes.csv")

write.csv(amastM_markers,
          "/datastore/Ross/cruzi_paper/bulk/paper/DEG_files/stageFine/amastMiddle_unique_sigGenes.csv")

write.csv(amastL_markers,
          "/datastore/Ross/cruzi_paper/bulk/paper/DEG_files/stageFine/amastLate_unique_sigGenes.csv")




amastE_df <- as.data.frame(amastE_markers)
amastE_logfc_order <- amastE_df[order(amastE_df$log2FoldChange, decreasing = T),]

amastM_df <- as.data.frame(amastM_markers)
amastM_logfc_order <- amastM_df[order(amastM_df$log2FoldChange, decreasing = T),]

amastL_df <- as.data.frame(amastL_markers)
amastL_logfc_order <- amastL_df[order(amastL_df$log2FoldChange, decreasing = T),]

epimast_df <- as.data.frame(epimast_markers)
epimast_logfc_order <- epimast_df[order(epimast_df$log2FoldChange, decreasing = T),]

trypomast_df <- as.data.frame(trypomast_markers)
trypomast_logfc_order <- trypomast_df[order(trypomast_df$log2FoldChange, decreasing = T),]

meta_df <- as.data.frame(metacyclic_markers)
meta_logfc_order <- meta_df[order(meta_df$log2FoldChange, decreasing = T),]

dim(amastE_df)

#Heatmap
sigGenes <- unique(c(row.names(amastE_logfc_order)[1:5], row.names(amastM_logfc_order)[1:5] ,
                     row.names(amastL_logfc_order)[1:5], row.names(trypomast_logfc_order)[1:5],  
                     row.names(epimast_logfc_order)[1:5], row.names(meta_logfc_order)[1:5]))


ntd <- normTransform(dds_cruzi)
library("pheatmap")
df <- as.data.frame(colData(dds_cruzi)[,c("stage_fine","stage_coarse")])

mtx_test <- as.matrix(assay(ntd)[sigGenes,])

row.subsections <- c(5,5,5,5, 5, 5)
row_split = data.frame("order"= rep(c("Early Amastigote markers", "Middle Amastigote markers", "Late Amastigote markers", "Trypomastigote markers","Epimastigote markers", "Metacyclic markers"), row.subsections))
row_split$order <- factor(row_split$order, c("Early Amastigote markers", "Middle Amastigote markers", "Late Amastigote markers", "Trypomastigote markers","Epimastigote markers", "Metacyclic markers"))

ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows = F,
                        row_split = row_split, cluster_row_slices = F)

mtx_test <- t(scale(t(as.matrix(assay(ntd)[sigGenes,]))))

pdf("/datastore/Ross/cruzi_paper/bulk/paper/plots/BulkRNA_deseq2_top10_fine_lifecycleMarkers_scaledHeatmap.pdf",
    height = 10)
ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows = F,
                        row_split = row_split, cluster_row_slices = F,
                        heatmap_legend_param = list(title = "Scaled expression"))
dev.off()




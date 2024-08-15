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

#set.tempdir("pooling-scRNA/")
ulimit::memory_limit(50000)

setwd("/datastore/Ross")

de_dir <- "/datastore/Ross/cruzi_files/lifecycle/bulk_RNA_seq/output/DE_tables/"

# gtf <- rtracklayer::import("cruzi_files/lifecycle/Tcruzi_Dm28c_2500UTR_MT_Hsampien_transcriptome.gtf")
# gtf_df=as.data.frame(gtf)



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
                                    design = ~0 + stage_fine)

dds_cruzi <- collapseReplicates(dds_cruzi, dds_cruzi$sample)

smallestGroupSize <- 11
keep <- rowSums(counts(dds_cruzi) >= 10) >= smallestGroupSize
dds_cruzi <- dds_cruzi[keep,]

dds_cruzi <- DESeq(dds_cruzi)
res <- results(dds_cruzi)






#Plot cell cycle genes
g1_markers <- read.delim("/datastore/Ross/cruzi_paper/bulk/input_objects/chavez_2017_data/dm28cOrthologs_of_chavez_cruzi_cellCycle_genes_syntenic/G1.csv", sep = ",")
g2_markers <- read.delim("/datastore/Ross/cruzi_paper/bulk/input_objects/chavez_2017_data/dm28cOrthologs_of_chavez_cruzi_cellCycle_genes_syntenic/G2.csv", sep=",")
s_markers <- read.delim("/datastore/Ross/cruzi_paper/bulk/input_objects/chavez_2017_data/dm28cOrthologs_of_chavez_cruzi_cellCycle_genes_syntenic/S.csv", sep=",")

#combine all markers into one list
cellcyclemarkers <- c(g1_markers$Gene.ID, g2_markers$Gene.ID, s_markers$Gene.ID) 
head(cellcyclemarkers)

cellcyclemarkers <- intersect(cellcyclemarkers, row.names(dds_cruzi))

g1_markers <- subset(g1_markers, Gene.ID %in% cellcyclemarkers)
g2_markers <- subset(g2_markers, Gene.ID %in% cellcyclemarkers)
s_markers <- subset(s_markers, Gene.ID %in% cellcyclemarkers)

#combine all markers into one list
#Taking only 50 genes due to how many we have to plot
cellcyclemarkers <- c(g1_markers$Gene.ID[1:40], g2_markers$Gene.ID[1:40], s_markers$Gene.ID[1:40]) 
head(cellcyclemarkers)

cellcyclemarkers <- intersect(cellcyclemarkers, row.names(dds_cruzi))


library("pheatmap")
ntd <- normTransform(dds_cruzi)
df <- as.data.frame(colData(dds_cruzi)[,c("stage_fine", "sample")])


df_test <- melt(assay(ntd)[cellcyclemarkers,])
mtx_test <- as.matrix(assay(ntd)[cellcyclemarkers,])

library(ComplexHeatmap)

row.subsections <- c(40,40,40)
row_split = data.frame(rep(c("G1", "G2M", "S"), row.subsections))

ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows = T,
                        row_split = row_split, cluster_row_slices = F)

mtx_test <- t(scale(t(mtx_test)))

pdf("/datastore/Ross/cruzi_paper/bulk/paper/plots/chavez_subset_cc_lifecycle_stages_heatmap.pdf", height = 22)
ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows =T,
                        row_split = row_split, cluster_row_slices = F,
                        heatmap_legend_param = list(title = "Scaled expression"))
dev.off()



#Plot cell cycle genes - Mapped not orthologged
g1_markers <- read.delim("/datastore/Ross/cruzi_paper/bulk/input_objects/chavez_2017_data/dm28cMapped_of_chavez_cruzi_cellCycle_genes/g1_genes_chavez_dm28c.csv", sep = ",")
g2_markers <- read.delim("/datastore/Ross/cruzi_paper/bulk/input_objects/chavez_2017_data/dm28cMapped_of_chavez_cruzi_cellCycle_genes/g2m_genes_chavez_dm28c.csv", sep=",")
s_markers <- read.delim("/datastore/Ross/cruzi_paper/bulk/input_objects/chavez_2017_data/dm28cMapped_of_chavez_cruzi_cellCycle_genes/s_genes_chavez_dm28c.csv", sep=",")

#combine all markers into one list
cellcyclemarkers <- c(g1_markers$x, g2_markers$x, s_markers$x) 
head(cellcyclemarkers)

cellcyclemarkers <- intersect(cellcyclemarkers, row.names(dds_cruzi))

g1_markers <- subset(g1_markers, x %in% cellcyclemarkers)
g2_markers <- subset(g2_markers, x %in% cellcyclemarkers)
s_markers <- subset(s_markers, x %in% cellcyclemarkers)

#combine all markers into one list
#Taking only 50 genes due to how many we have to plot
cellcyclemarkers <- c(g1_markers$x[1:98], g2_markers$x[1:73], s_markers$x[1:80]) 
head(cellcyclemarkers)

cellcyclemarkers <- intersect(cellcyclemarkers, row.names(dds_cruzi))


library("pheatmap")
ntd <- normTransform(dds_cruzi)
df <- as.data.frame(colData(dds_cruzi)[,c("stage_fine", "sample")])


df_test <- melt(assay(ntd)[cellcyclemarkers,])
mtx_test <- as.matrix(assay(ntd)[cellcyclemarkers,])

library(ComplexHeatmap)

row.subsections <- c(98,73,80)
row_split = data.frame(rep(c("G1", "G2M", "S"), row.subsections))

ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows = T,
                        row_split = row_split, cluster_row_slices = F)

mtx_test <- t(scale(t(mtx_test)))

pdf("/datastore/Ross/cruzi_paper/bulk/paper/plots/chavezMapped_subset_cc_lifecycle_stages_heatmap.pdf", height = 10)
ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows =T,
                        row_split = row_split, cluster_row_slices = F,
                        heatmap_legend_param = list(title = "Scaled expression"),
                        row_labels = rep("", dim(mtx_test)[1]))
dev.off()




















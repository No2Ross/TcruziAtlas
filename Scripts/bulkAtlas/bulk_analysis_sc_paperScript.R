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




#Plot single cell markers (all cell types)
ntd <- normTransform(dds_cruzi)

sc_markers <- read.csv("cruzi_files/lifecycle/clean_analysis/top_markers_clusters_no48hr_filtered.csv")
sc_markers <- sc_markers[,-1]
sc_markers <- sc_markers[order(sc_markers$logfoldchanges, decreasing = T),]

marker_count <- c()
marker_genes <- c()

group_list <- c("early_mid_amast_5", "late_amast_1", "trypomast_0", "trypomast_4", "trypomast_9", "epimast_2",
                "epi_meta_trans_6", "epi_meta_trans_7", "meta_3", "meta_8")

for(i in group_list){
  current <- subset(sc_markers, group == i)
  
  current <- current$names
  current <- intersect(current, row.names(dds_cruzi))
  
  if(length(current) < 5){
    marker_count <- append(marker_count, length(current))
    marker_genes <- append(marker_genes, current)
  }
  
  else{
    marker_count <- append(marker_count, 5)
    marker_genes <- append(marker_genes, current[1:5])
  }
  
  
}

mtx_test <- as.matrix(assay(ntd)[marker_genes,])

row.subsections <- marker_count
row_split = data.frame("cluster" = rep( paste(group_list, "\n markers"), row.subsections))
row_split$cluster <- factor(row_split$cluster, paste(group_list, "\n markers"))

ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows = F,
                        row_split = row_split, cluster_row_slices = F,
                        row_title_rot = switch("left", "left" = 0, "right" = 270))

mtx_test <- t(scale(t(mtx_test)))

pdf("/datastore/Ross/cruzi_paper/bulk/paper/plots/BulkRNAseq_lifecycle_stages_scMarkers_heatmap.pdf", 
    height = 10)
ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows=F,
                        row_split = row_split, cluster_row_slices = F,
                        row_title_rot = switch("left", "left" = 0, "right" = 270),
                        heatmap_legend_param = list(title = "Scaled expression"))
dev.off()



#Metacyclogenesis
epi_meta_markers <- read.csv("/datastore/Ross/cruzi_paper/single_cell/epi_meta_trans_DE_genes_peak.csv", row.names = 1)
epi_meta_markers <- epi_meta_markers[order(epi_meta_markers$meanLogFC, decreasing = T),]

row.names(epi_meta_markers) <- str_replace_all(row.names(epi_meta_markers), "-", "_")

early_markers <- subset(epi_meta_markers, cluster == "early")
early_markers <- row.names(early_markers[-which(early_markers %in% c("MT tcruzi", "MT2tcruzi"))])
early_markers <- early_markers[1:10]

middle_markers <- subset(epi_meta_markers, cluster == "middle")
middle_markers <- row.names(middle_markers[-which(middle_markers %in% c("MT tcruzi", "MT2tcruzi"))])
middle_markers <- middle_markers[1:10]

late_markers <- subset(epi_meta_markers, cluster == "late")
late_markers <- row.names(late_markers[-which(late_markers %in% c("MT tcruzi", "MT2tcruzi"))])
late_markers <- late_markers[1:10]

sigGenes_epi_meta <- c(early_markers, middle_markers, late_markers)

epi_meta_dds <- dds_cruzi[, -which(colnames(dds_cruzi) %in% c("EarlyAmast1","EarlyAmast2", "MidAmast1","MidAmast2", 
                                                                            "LateAmast1","LateAmast2","Trypomast1", "Trypomast2"))]

ntd <- normTransform(epi_meta_dds)

sigGenes_epi_meta <- intersect(sigGenes_epi_meta, row.names(ntd))

mtx_test <- as.matrix(assay(ntd)[sigGenes_epi_meta,])

row.subsections <- c(10,10,10)
row_split = data.frame("order"= rep(paste(c("Early", "Middle", "Late"), "\n markers"), row.subsections))
row_split$order <- factor(row_split$order, paste(c("Early", "Middle", "Late"), "\n markers"))

ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows = F,
                        row_split = row_split)

mtx_test <- t(scale(t(as.matrix(assay(ntd)[sigGenes_epi_meta,]))))

pdf("/datastore/Ross/cruzi_paper/bulk/paper/plots/epi_to_meta_DEGs_scaled_heatmap.pdf", height = 7)
ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows = F,
                        row_split = row_split, cluster_row_slices = F,
                        heatmap_legend_param = list(title = "Scaled expression"))
dev.off()








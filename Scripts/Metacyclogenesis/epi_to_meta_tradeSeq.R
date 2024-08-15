#This analysis doesn't exclude the meta_8 cells. Remove them and redo analysis


library(pheatmap)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
library(metR)


library(unixtools)

#set.tempdir("pooling-scRNA/")
ulimit::memory_limit(60000)

setwd("/datastore/Ross")

raw_mtx <- read.table("/datastore/Ross/cruzi_paper/single_cell/input_objects/epi_meta/epi_to_meta_expMtxRaw.csv", sep = ",")
genes_frame <- read.csv("/datastore/Ross/cruzi_paper/single_cell/input_objects/epi_meta/epi_to_meta_genes.csv", row.names = 1)
metaData <-  read.csv("/datastore/Ross/cruzi_paper/single_cell/input_objects/epi_meta/epi_to_meta_metadata.csv")

raw_mtx <- t(raw_mtx)

colnames(raw_mtx) <- metaData$X
row.names(raw_mtx) <- genes_frame[,1]
row.names(metaData) <- metaData$X

# norm_mtx <- read.table("/datastore/Ross/cruzi_paper/single_cell/input_objects/epi_meta/epi_to_meta_expMtx.csv", sep = ",")
# norm_mtx <- t(norm_mtx)
# row.names(norm_mtx) <- genes_frame[,1]
# colnames(norm_mtx) <- metaData$X


epi_meta <- CreateSeuratObject(raw_mtx, meta.data = metaData)

#We are interested in the transition between epi and meta
#Because of this we want to get rid of cells at the very metacyclic end of things
#Only keep cells with pseudotime < 0.6

# epi_meta@assays$RNA@data <- as.sparse(norm_mtx)

hist(epi_meta$dpt_pseudotime)
epi_meta <- subset(epi_meta, dpt_pseudotime < 0.65)
# epi_meta <- subset(epi_meta, stages_subset != "meta8")

epi_meta_sce <- as.SingleCellExperiment(epi_meta)
rm(raw_mtx)
rm(metaData)
gc()


#tradeseq
# set.seed(3)
# icMat <- evaluateK(counts = as.matrix(assays(epi_meta_sce)$counts),
#                    pseudotime = epi_meta_sce$dpt_pseudotime,
#                    cellWeights = rep(1, length(epi_meta_sce$dpt_pseudotime)),
#                    nGenes = 300,
#                    k = 4:12)
# 
# rm(icMat)
# gc()

cruzi_GAM <- fitGAM(as.matrix(assays(epi_meta_sce)$counts), pseudotime = epi_meta_sce$dpt_pseudotime, cellWeights = rep(1, length(epi_meta_sce$dpt_pseudotime)), nknots = 7, parallel = F)

cruzi_GAM$slingshot$stage_cluster <- epi_meta_sce$stages_subset

saveRDS(cruzi_GAM, "/datastore/Ross/cruzi_paper/single_cell/output_objects/epi_meta_trans_tradeSeq.rds")
cruzi_GAM <- readRDS("/datastore/Ross/cruzi_paper/single_cell/output_objects/epi_meta_trans_tradeSeq.rds")

library(stringr)

assoRes <- associationTest(cruzi_GAM, l2fc = 0.25)

assoRes$p_adjust <- p.adjust(assoRes$pvalue, method = "bonferroni")
assoRes_sig <- subset(assoRes, p_adjust < 0.05)

#Remove genes which are expressed in less than 2.5% of the cells
gene_stats <- rowSums(as.matrix(cruzi_GAM@assays@data@listData$counts) > 0)
genes_keep <- names(gene_stats)[which(gene_stats > ( dim(cruzi_GAM)[1] * 0.1 ) )]

assoRes_sig_filter <- assoRes_sig[intersect(genes_keep,row.names(assoRes_sig)),]

assoRes_sig_filter <- assoRes_sig_filter[str_detect(row.names(assoRes_sig_filter), "C4B63"),]

write.csv(assoRes_sig_filter, "/datastore/Ross/cruzi_paper/single_cell/DEG_files/tradeSeq_epi_meta_filtered_DEgenes.csv")
write.csv(assoRes, "/datastore/Ross/cruzi_paper/single_cell/DEG_files/tradeSeq_epi_meta_nonFiltered_DEgenes.csv")

yhatSmooth <- predictSmooth(cruzi_GAM, gene = row.names(assoRes_sig_filter), nPoints = 50, tidy = FALSE)

x <- t(scale(t(log1p(cruzi_GAM@assays@data$counts[row.names(assoRes_sig_filter),]))))

# test <- pheatmap(x,
#                        cluster_cols = FALSE,
#                        show_rownames = FALSE, 
#                        show_colnames = FALSE,
#                        clustering_distance_rows = "correlation")

heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                       cluster_cols = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE,
                       clustering_distance_rows = "correlation")

cl <- cutree(heatSmooth$tree_row, k = 4)

gene_groups <- data.frame("Cluster" = cl, 
                          row.names = names(cl))

#andrew wave heatmap
source("/datastore/Ross/cruzi_paper/single_cell/scripts/andrew_waves.R")
sce <- SingleCellExperiment( list( counts = as.sparse(t(scale(t(yhatSmooth[, 1:50])))), data = as.sparse(t(scale(t(yhatSmooth[, 1:50])))) ) )
sce$pseudotime <- seq(1, 50, 1)
waves <- get_waves(sce, "pseudotime")

waves <- waves[order(waves$phase),]

reorder_test <- yhatSmooth[order(match(row.names(yhatSmooth),row.names(waves))),]


#Cluster the genes according to their phase: 0-16 = early, 17 - 33 = middle, 34-50 = late
mtx_test <- t(scale(t(reorder_test[, 1:50])))

gene_cluster <- data.frame("cluster" = rep("na", dim(waves)[1]) , "phase" = waves$phase, row.names = waves$gene)

counter <- 0
for(i in 1:length(gene_cluster$phase)){
  
  current_phase <- gene_cluster$phase[i]
  
  if( current_phase <= 16){
    gene_cluster$cluster[i] <- "early"
  }
  
  else if( current_phase > 16 & current_phase <=29){
    counter <- counter + 1
    gene_cluster$cluster[i] <- "middle"
    
  }
  
  else{
    gene_cluster$cluster[i] <- "late"
  }
  
  
}

cluster_label_number <- table(gene_cluster$cluster)

gene_cluster <- gene_cluster[order(match(row.names(gene_cluster),row.names(assoRes_sig_filter))),]

gene_cluster$meanLogFC <- assoRes_sig_filter$meanLogFC
gene_cluster$pvalue <- assoRes_sig_filter$pvalue
gene_cluster$p_adjust <- assoRes_sig_filter$p_adjust

write.csv(gene_cluster, "/datastore/Ross/cruzi_paper/single_cell/epi_meta_trans_DE_genes_peak.csv")
gene_cluster <- read.delim("/datastore/Ross/cruzi_paper/single_cell/epi_meta_trans_DE_genes_peak.csv", sep = ",", header = T)

library(ComplexHeatmap)

row.subsections <- c(cluster_label_number["early"],cluster_label_number["middle"],cluster_label_number["late"])
row_split = data.frame("order"= rep(c("Early peak", "Middle peak", "Late peak"), row.subsections))
row_split$order <- factor(row_split$order, c("Early peak", "Middle peak", "Late peak"))

ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows = F,
                        row_split = row_split, cluster_row_slices = F,
                        heatmap_legend_param = list(title = "Scaled Smoothed expression"),
                        row_labels = rep("", dim(mtx_test)[1]), column_labels = rep("", dim(mtx_test)[2]),
                        column_title = "Pseudotime", row_title_rot = 0)

pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/tradeSeq_genes_heatmap.pdf")
ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows = F,
                        row_split = row_split, cluster_row_slices = F,
                        heatmap_legend_param = list(title = "Scaled expression"),
                        row_labels = rep("", dim(mtx_test)[1]), column_labels = rep("", dim(mtx_test)[2]),
                         column_title = "Pseudotime", row_names_rot = 180)
dev.off()

assoRes_sig_filter_phase <- assoRes_sig_filter
assoRes_sig_filter_phase <- assoRes_sig_filter_phase[order(match(row.names(assoRes_sig_filter_phase), row.names(gene_cluster))),]
assoRes_sig_filter_phase$phase_cluster <- gene_cluster$cluster

#Early peak genes
#2-amino-3-ketobutyrate coenzyme A ligase
i <- "C4B63-32g214"
j <- "2-amino-3-ketobutyrate coenzyme A ligase"
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/early_peak/tradeSeq_genes_C4B63_32g214.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Threonine-dehydrogenase
i <- "C4B63-42g129"
j <- "L-Threonine-dehydrogenase"
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/early_peak/tradeSeq_genes_C4B63_12g175.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


#Proline racemase
i <- "C4B63-52g138"
j <- "Proline racemase"
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/early_peak/tradeSeq_genes_C4B63_52g138.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Trans-sialidase group II
#This paper suggests that in culture stationary epimastigotes express trans-sialidase- Trans-sialidase from Trypanosoma cruzi epimastigotes is expressed at the stationary phase and 
#is different from the enzyme expressed in trypomastigotes

#(Trans-sialidase from Trypanosoma cruzi epimastigotes is expressed at the stationary phase and is different from the enzyme expressed in trypomastigotes)

i <- "C4B63-209g8"
j <- "Trans-sialidase group II"
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/early_peak/tradeSeq_genes_C4B63_209g8.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

i <- "C4B63-7g438"
j <- "Glutamate Dehydrogenase"
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/early_peak/tradeSeq_genes_C4B63_7g438.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


#middle peak

i <- "C4B63-12g175"
j <- "Protein Associated with Differentiation"
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/mid_peak/tradeSeq_genes_C4B63_12g175.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

i <- "C4B63-83g24"
j <- "ZC3H12"
#A Trypanosoma cruzi zinc finger protein that is implicated in the control of epimastigote-specific gene expression and metacyclogenesis
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/mid_peak/tradeSeq_genes_C4B63_83g24.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Might not be PBP1, but a PBP1 binding protein
i <- "C4B63-21g42"
j <- "PBP1"
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/mid_peak/tradeSeq_genes_C4B63_21g42.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

i <- "C4B63-18g156"
j <- "MKT1"
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/mid_peak/tradeSeq_genes_C4B63_18g156.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

i <- "C4B63-34g1526c"
j <- "TcSMUG-L"
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/mid_peak/tradeSeq_genes_C4B63_34g1526c.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

i <- "C4B63-6g507"
j <- "RAD51"
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/mid_peak/tradeSeq_genes_C4B63_6g507.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()




#Late peak

#Citrate synthase
i <- "C4B63-4g184"
j <- "Citrate Synthase"
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/late_peak/tradeSeq_genes_C4B63_4g184.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


i <- "C4B63-16g183"
j <- "surface protease GP63"
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/late_peak/tradeSeq_genes_C4B63_16g183.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


i <- "C4B63-18g326"
j <- "ZFP1"
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/late_peak/tradeSeq_genes_C4B63_18g326.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

i <- "C4B63-32g1409c"
j <- "Chagasin"
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/late_peak/tradeSeq_genes_C4B63_32g1409c.pdf", width = 9)
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = i,
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle(paste0(i," \n (", j, ")")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()



ggplot() + geom_point(aes(epi_meta_sce$dpt_pseudotime , as.matrix(assays(epi_meta_sce)$counts)["C4B63-35g186",]))
ggplot() + geom_point(aes(colSums(as.matrix(assays(epi_meta_sce)$counts)), epi_meta_sce$dpt_pseudotime))

#Epimast
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-81g59", pointCol = cruzi_GAM$slingshot$stage_cluster)+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))




#RNA Binding genes

#Genes which are not commented on in Cruzi literature
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_new_C4B63_27g1020c.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-27g1020c",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_27g1020c") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_new_C4B63_58g143.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-58g143",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_58g143") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


#Marta list
#MT	Minning, 2009; Garcia-Huertas, 2022
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_C4B63-162g31_MT.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-162g31",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_162g31") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#MT	Minning, 2009; Garcia-Huertas, 2022
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_C4B63-46g10_MT.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-46g10",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_46g10") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#EP	Tavares, 2020
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_C4B63-9g319_EP.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-9g319",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_9g319") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#TCT, MT	Li, 2016; Tavares, 2020; Belew, 2017; Sabalette, 2023; Santos, 2018
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_C4B63-23g237_TCT_MT.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-23g237",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_23g237") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Metacyclogenesis	Sabalette, 2013
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_C4B63-19g103_MT.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-19g103",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_19g103") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#AMA, EP	Li, 2016; Tavares, 2020
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_C4B63-7g426_AMA_EP.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-7g426",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_7g426") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#MT and TCT	Li, 2016; Sabalette, 2023; Garcia-Huertas, 2022; Tavares, 2020
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_C4B63-27g167_MT_TCT.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-27g167",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_27g167") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#MT	Minning, 2009
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_C4B63-4g545_MT.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-4g545",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_4g545") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#MT and early AMA	Li, 2016; Sabalette, 2023; Santos, 2018
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_C4B63-16g295_MT_earlyAMA.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-16g295",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_16g295") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#EP	Tavares, 2020; Santos, 2018
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_C4B63-18g216_EP.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-18g216",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_18g216") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#EP	Tavares, 2020
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_C4B63-2g318_EP.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-2g318",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_2g318") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#EP	Tavares, 2020

pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_C4B63-6g149_EP.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-6g149",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_6g149") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#EP	Tavares, 2020
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_C4B63-34g364_EP.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-34g364",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_34g364") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#AMA, EP	Tavares, 2020; Li, 2016
pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_C4B63-6g592_EP_AMA.pdf")
plotSmoothers(cruzi_GAM, as.matrix(assays(epi_meta_sce)$counts) , gene = "C4B63-6g592",
              pointCol = cruzi_GAM$slingshot$stage_cluster) + ggtitle("C4B63_6g592") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

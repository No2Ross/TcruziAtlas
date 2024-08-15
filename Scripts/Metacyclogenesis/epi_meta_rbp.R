library(pheatmap)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
library(metR)
library(stringr)

library(unixtools)

#set.tempdir("pooling-scRNA/")
ulimit::memory_limit(60000)

setwd("/datastore/Ross")

cruzi_GAM <- readRDS("/datastore/Ross/cruzi_paper/single_cell/output_objects/epi_meta_trans_tradeSeq.rds")

rna_binding_genes <- c(
  "C4B63_84g49", "C4B63_23g237", "C4B63_113g39", "C4B63_31g222", "C4B63_10g157",
  "C4B63_32g47", "C4B63_162g31", "C4B63_46g10", "C4B63_45g126", "C4B63_21g94",
  "C4B63_19g123", "C4B63_80g35", "C4B63_80g36", "C4B63_46g24", "C4B63_27g1020c",
  "C4B63_42g189", "C4B63_20g198", "C4B63_9g319", "C4B63_50g198", "C4B63_8g299",
  "C4B63_13g332", "C4B63_10g1625c", "C4B63_44g179", "C4B63_17g180", "C4B63_13g1743c",
  "C4B63_13g331", "C4B63_13g1746c", "C4B63_72g542c", "C4B63_38g46", "C4B63_30g219",
  "C4B63_23g237", "C4B63_52g148", "C4B63_32g272", "C4B63_18g205", "C4B63_19g191",
  "C4B63_10g449", "C4B63_19g103", "C4B63_19g101", "C4B63_19g98", "C4B63_8g490",
  "C4B63_8g492", "C4B63_63g6", "C4B63_7g426", "C4B63_4g356", "C4B63_10g483",
  "C4B63_2g50", "C4B63_2g154", "C4B63_2g155", "C4B63_7g223", "C4B63_55g172",
  "C4B63_4g172", "C4B63_58g143", "C4B63_27g167", "C4B63_8g41", "C4B63_4g545",
  "C4B63_21g53", "C4B63_42g82", "C4B63_9g235", "C4B63_16g295", "C4B63_57g971c",
  "C4B63_18g278", "C4B63_17g65", "C4B63_17g64", "C4B63_297g9", "C4B63_17g54",
  "C4B63_17g52", "C4B63_12g28", "C4B63_23g29", "C4B63_2g96", "C4B63_20g41",
  "C4B63_7g282", "C4B63_45g133", "C4B63_3g1103", "C4B63_2g768", "C4B63_18g216",
  "C4B63_2g331", "C4B63_2g174", "C4B63_17g244", "C4B63_76g21", "C4B63_206g13",
  "C4B63_13g223", "C4B63_21g102", "C4B63_20g203", "C4B63_69g28", "C4B63_15g5",
  "C4B63_397g11", "C4B63_13g208", "C4B63_2g567", "C4B63_4g429", "C4B63_4g427",
  "C4B63_12g213", "C4B63_95g32", "C4B63_83g24", "C4B63_2g33", "C4B63_2g232",
  "C4B63_76g31", "C4B63_76g18", "C4B63_76g17", "C4B63_2g343", "C4B63_2g319",
  "C4B63_2g318", "C4B63_6g149", "C4B63_6g105", "C4B63_18g169", "C4B63_54g10",
  "C4B63_128g60", "C4B63_116g22", "C4B63_12g400", "C4B63_97g13", "C4B63_97g14",
  "C4B63_97g15", "C4B63_8g463", "C4B63_112g44", "C4B63_252g14", "C4B63_2g507",
  "C4B63_10g480", "C4B63_10g480", "C4B63_7g312", "C4B63_133g32", "C4B63_97g1",
  "C4B63_4g276", "C4B63_4g275", "C4B63_12g2035c", "C4B63_12g327", "C4B63_34g364",
  "C4B63_14g59", "C4B63_30g164", "C4B63_122g23", "C4B63_27g6", "C4B63_18g330",
  "C4B63_18g326", "C4B63_18g253", "C4B63_18g256", "C4B63_142g4", "C4B63_114g38",
  "C4B63_238g11", "C4B63_12g395", "C4B63_63g14", "C4B63_4g1393c", "C4B63_34g360",
  "C4B63_10g457", "C4B63_42g230", "C4B63_78g39", "C4B63_2g787", "C4B63_7g407",
  "C4B63_13g170", "C4B63_43g137", "C4B63_40g132", "C4B63_22g175", "C4B63_70g111",
  "C4B63_70g112", "C4B63_11g490", "C4B63_13g295", "C4B63_6g592", "C4B63_97g23",
  "C4B63_17g243", "C4B63_41g148", "C4B63_413g14", "C4B63_439g12", "C4B63_13g170",
  "C4B63_13g295", "C4B63_58g56", "C4B63_25g296", "C4B63_25g296", "C4B63_25g297",
  "C4B63_72g94", "C4B63_208g34", "C4B63_208g34", "C4B63_92g63", "C4B63_232g18",
  "C4B63_14g144", "C4B63_18g225", "C4B63_51g844c", "C4B63_51g842c",
  "C4B63_51g842c", "C4B63_7g209", "C4B63_7g210", "C4B63_42g195", "C4B63_35g189",
  "C4B63_58g123", "C4B63_10g153", "C4B63_10g154", "C4B63_112g27", "C4B63_11g39",
  "C4B63_13g318", "C4B63_148g34", "C4B63_153g280c", "C4B63_187g2", "C4B63_189g33",
  "C4B63_18g331", "C4B63_18g335", "C4B63_204g19", "C4B63_227g8", "C4B63_278g12",
  "C4B63_278g8", "C4B63_28g171", "C4B63_293g10", "C4B63_40g145", "C4B63_40g185",
  "C4B63_49g63", "C4B63_4g172", "C4B63_4g382", "C4B63_4g551", "C4B63_56g51",
  "C4B63_8g515"
)

rna_binding_genes <- str_replace_all(rna_binding_genes, "_", "-")

rna_binding_genes <- intersect(rna_binding_genes, row.names(cruzi_GAM))

gene_cluster <- read.delim("/datastore/Ross/cruzi_paper/single_cell/epi_meta_trans_DE_genes_peak.csv", sep = ",", row.names = 1)

row.names(gene_cluster) <- str_replace_all(row.names(gene_cluster), "_", "-")

rna_binding_genes <- intersect(rna_binding_genes, row.names(gene_cluster))

yhatSmooth <- predictSmooth(cruzi_GAM, gene = rna_binding_genes, nPoints = 50, tidy = FALSE)

rna_binding_df <- gene_cluster[rna_binding_genes,]

rna_binding_df_save <- rna_binding_df

row.names(rna_binding_df_save) <- str_replace_all(row.names(rna_binding_df_save), "-", "_")

write.csv(rna_binding_df_save, "/datastore/Ross/cruzi_paper/single_cell/rna_binding_genes_DE.csv")

#Remove genes which are not expressed in at least 10% of the cells
count_mtx <- 1 - (rowCounts(cruzi_GAM@assays@data@listData$counts[rna_binding_genes,], value = 0, useNames = T) /length(colnames(cruzi_GAM@assays@data@listData$counts[rna_binding_genes,])))

valid_rna_binding_genes <- names(count_mtx)[which(count_mtx > 0.1)]

yhatSmooth <- predictSmooth(cruzi_GAM, gene = valid_rna_binding_genes, nPoints = 50, tidy = FALSE)

heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                       cluster_rows = F,
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

cluster_label_number <- table(rna_binding_df$cluster)

#Cluster the genes according to their phase: 0-16 = early, 17 - 33 = middle, 34-50 = late
mtx_test <- t(scale(t(reorder_test[, 1:50])))

row.subsections <- c(cluster_label_number["early"],cluster_label_number["middle"],cluster_label_number["late"])
row_split = data.frame("order"= rep(c("Early peak RBPs", "Middle peak RBPs", "Late peak RBPs"), row.subsections))
row_split$order <- factor(row_split$order, c("Early peak RBPs", "Middle peak RBPs", "Late peak RBPs"))


pdf("/datastore/Ross/cruzi_paper/single_cell/plots/epi_meta/RNA_binding_heatmap.pdf")
ComplexHeatmap::Heatmap(mtx_test, cluster_columns = F, cluster_rows = F,
                        row_split = row_split, cluster_row_slices = F,
                        row_labels =str_replace_all(row.names(mtx_test), "-", "_") ,
                        heatmap_legend_param = list(title = "Scaled Smoothed expression"), column_labels = rep("", dim(mtx_test)[2]),
                        column_title = "Pseudotime", row_title_rot = 0)
dev.off()







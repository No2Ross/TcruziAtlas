get_waves <- function(
  sce,
  pseudotime_slot = "slingPseudotime_1",
  n_cores = 1
) {
  
  if ( !any(colnames(sce@colData) == pseudotime_slot)) {
    stop(paste0("Pseudotime slot '", pseudotime_slot ,"' does not exist"))
  }
  pseudotime = sce@colData[[pseudotime_slot]]
  
  heatmap_counts = SingleCellExperiment::counts(sce)[,order(pseudotime)]
  
  # TODO AM rewrite with https://bioconductor.org/packages/release/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.html#single-machine
  waves_list = parallel::mclapply(rownames(heatmap_counts), function (gene) {
    wave = as.data.frame(FitWave(as.matrix(heatmap_counts[gene,]), 1))
    rownames(wave) = c(gene)
    return(wave)
  }, mc.cores=n_cores)
  waves <- do.call("rbind", waves_list)
  
  waves = as.data.frame(waves)
  waves$phase = waves$phase * (((180/3.141593)/360)*max(pseudotime)) # pseudotime
  waves$gene = rownames(waves)
  
  # Add power to waves
  waves$total_expression = rowSums(heatmap_counts[rownames(waves),])
  # Bozdech et al. 2003 (https://doi.org/10.1371/journal.pbio.0000005) use plus
  # or minus 1/48 (i.e. one bulk either side of the peak, 6.25% window)
  # Here we use plus or minus 5% (i.e. a 10% window)
  five_percent_of_pdt = 0.05*max(waves$phase)
  waves$peak_expression = 0
  for (gene in rownames(waves)) {
    waves[gene,"peak_expression"] = sum(heatmap_counts[gene,pseudotime > waves[gene,]$phase-five_percent_of_pdt & pseudotime < waves[gene,]$phase+five_percent_of_pdt])
    waves[gene,"cellcount_in_peak"] = length(heatmap_counts[gene,pseudotime > waves[gene,]$phase-five_percent_of_pdt & pseudotime < waves[gene,]$phase+five_percent_of_pdt])
  }
  waves$power = waves$amplitude / waves$peak_expression
  
  return(waves)
}


gene_selection_matrix <- function(x, waves, genes=c(), pseudotime_slot="slingPseudotime_1", matrix_rows=1000, matrix_cols=1000, n_cores=1){
  
# gene_selection_matrix <- function(x, waves, genes=c(), pseudotime_slot="slingPseudotime_1", target_matrix_size=1000, n_cores=1){
  
  if ( !any(colnames(x@colData) == pseudotime_slot)) {
    stop(paste0("Pseudotime slot '", pseudotime_slot ,"' does not exist"))
  }
  
  # First we need to subset only the requested genes
  if (length(genes) > 0) {
    # R passes parameters by value not reference so this is safe
    x = x[rownames(x) %in% genes,]
    waves = waves[rownames(waves) %in% genes,]
  }
  
  # Then get the normalised count matrix
  pseudotime = x@colData[[pseudotime_slot]]
  heatmap_counts = SingleCellExperiment::counts(x)[,order(pseudotime)]
  heatmap_counts = heatmap_counts[order(waves$phase),]
  
  small_heatmap_counts = redim_matrix(
    heatmap_counts,
    target_height = matrix_rows,
    target_width = matrix_cols,
    n_core=n_cores
  )
  
  heatmap_counts_ordered.df <- reshape2::melt(small_heatmap_counts, c("gene", "cell"), value.name = "log_expression")
  
  cell_sym = ggplot2::sym("cell")
  gene_sym = ggplot2::sym("gene")
  log_expr_sym = ggplot2::sym("log_expression")
  
  plot = ggplot2::ggplot(data=heatmap_counts_ordered.df,ggplot2::aes(x={{cell_sym}},y={{gene_sym}},fill={{log_expr_sym}})) +
    ggplot2::geom_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank()) +
    ggplot2::scale_fill_gradient(low = "white",
                                 high = "red",
                                 guide = "colorbar")
  
  return(plot)
}

redim_matrix <- function(
  mat,
  target_height = 100,
  target_width = 100,
  summary_func = function(x) mean(x, na.rm = TRUE),
  n_core = 1
) {
  
  if(target_height > nrow(mat) | target_width > ncol(mat)) {
    stop("Input matrix must be bigger than target width and height.")
  }
  
  seq_height <- round(seq(1, nrow(mat), length.out = target_height + 1))
  seq_width  <- round(seq(1, ncol(mat), length.out = target_width  + 1))
  
  # complicate way to write a double for loop
  # TODO AM rewrite with https://bioconductor.org/packages/release/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.html#single-machine
  do.call(rbind, parallel::mclapply(seq_len(target_height), function(i) { # i is row
    vapply(seq_len(target_width), function(j) { # j is column
      summary_func(
        mat[
          seq(seq_height[i], seq_height[i + 1]),
          seq(seq_width[j] , seq_width[j + 1] )
        ]
      )
    }, 0.0)
  }, mc.cores = n_core))
}

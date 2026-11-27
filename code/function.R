#' clean the rownames to only keep the gene symbol
clean_name <- function(x) sub(".*\\.", "", x)

#' Look up full rownames by gene symbol suffix
getrn <- function(se, x){
  grep(paste0("\\.", x, "$", collapse="|"), row.names(se), value=TRUE)
}

#' Calculate a "stress signature" coefficient for each sample
#'
#' For a given set of gene signature, quantify the degree of similarity 
#' between each sample's expression and the average response of stressed 
#' animals. 
#' This is done by regressing each sample's log2FC (over the signature genes) 
#' against the mean log2FC the Stress samples, using a robust (M-estimator) 
#' regression through the origin. The resulting slope is the "stress coefficient": 
#' ~1 means the sample looks like an average stressed animal, 
#' ~0 means no stress-like response.
#'
#' @param SE A SummarizedExperiment with a `log2FC` assay and colData
#'   columns `SampleName`, `Condition`, `Sex`, `group`.
#' @param gene_list Character vector of signature genes (rownames in `SE`).
#' @param dataset_name Label identifying the source of this calculation.
#'   The special value `"Reference\n"` triggers per-sample output (used
#'   for the reference/full-data analysis); any other value triggers
#'   per-group averaging (used when this function is called inside a
#'   permutation loop, since we only want one summary value per group
#'   per permutation, not one per animal).
#'
#' @return A data.frame. If `dataset_name == "Reference\n"`: one row per
#'   sample with columns group, Condition, Sex, coefficient, dataset.
#'   Otherwise: one row per group x Condition x Sex combination, with
#'   coefficient averaged across samples in that combination.
#' 
calculate_stress_coef <- function(SE, gene_list, dataset_name) {

  #For reference, calculate coef for all samples

  if(dataset_name == "Reference\n"){
    cols <- c("SampleName", "Condition", "Sex", "group")
    log2fc_mat <- assay(SE, "log2FC")
    e <- log2fc_mat[gene_list, ]
  
    stressed_samples <- SE$SampleName[SE$Condition == "Stress"]
    x <- rowMeans(e[, stressed_samples]) 
    # Here we calculate the mean log2FC of stressed animals as the reference. Every sample is later scored against this reference.
    
    # Calculate the rlm (robust regression) coefficient for each sample
    coefficients <- sapply(colnames(e), function(current_sample) {
      y <- e[, current_sample]
      
      # Model: y ~ 0 + x (regress sample's response against the average stressed response)
      model <- MASS::rlm(y ~ 0 + x)
      return(coef(model)[1])
    })
    df <- as.data.frame(colData(SE)) %>%
      dplyr::select(group, Condition, Sex) %>%
      mutate(
        coefficient = coefficients,
        dataset = dataset_name
      )
    return(df)

  } 

  # For baseline average the coef on a per group basis since permutations are computed later on

  else{
    cols <- c("SampleName", "Condition", "Sex", "group")
    log2fc_mat <- assay(SE, "log2FC")
    e <- log2fc_mat[gene_list, ]
  
    stressed_samples <- SE$SampleName[SE$Condition == "Stress"]
    x <- rowMeans(e[, stressed_samples])
    
    # Calculate the rlm coefficient for each sample
    coefficients <- sapply(colnames(e), function(current_sample) {
      y <- e[, current_sample]
      
      # Model: y ~ 0 + x (regress sample's response against the average stressed response)
      model <- MASS::rlm(y ~ 0 + x)
      return(coef(model)[1])
    })
  
    summary_df <- as.data.frame(colData(SE)) %>%
      dplyr::select(group, Condition, Sex) %>%
      mutate(
        coefficient = coefficients,
        dataset = dataset_name
      ) %>%
      group_by(dataset, group, Condition, Sex) %>%
      summarize(coefficient = mean(coefficient, na.rm = TRUE), .groups = "drop")
    return(summary_df)
  }
}

#' Plot sex-comparison views of a geneset's stress response
#'
#' Based on a given geneset (for example a set of stress DEGs), produces two
#' plots that compare how males and females respond to stress:
#'   1. A scatter plot of each gene's mean stress log2FC in males vs.
#'      females, to visualise whether the genes respond
#'      concordantly or divergently between the two sexes.
#'   2. A boxplot of each sample's "coefficient score": a score of how strongly 
#'      that sample's log2FC profile matches the average stressed-animal 
#'      profile, computed the same way as in `calculate_stress_coef()`), 
#'      split by group and Condition, with a Male vs Female t-test within each Condition.
#'
#' Note: gene names in `gene_list` and in `rownames(SE)` are expected to
#' be prefixed (e.g. "ENSMUSG00000012345.Hes5"); for the scatter plot
#' only, prefixes are stripped down to the gene symbol
#' (`sub(".*\\.", "", x)`) purely for cleaner point labels. The boxplot
#' coefficient calculation uses the original (unstripped) IDs, matching
#' how `gene_list` was generated upstream (e.g. from `stagedFDRwrapper`).
#'
#' @param SE A SummarizedExperiment with a `log2FC` assay and colData
#'   columns `SampleName`, `Sex`, `Condition`, `group`.
#' @param gene_list Character vector of the geneset (rownames in `SE`).
#' @param text1 Label used in both plot titles (e.g. "Interaction",
#'   describing which DEG set `gene_list` represents).
#' @param text2 Label used in the y-axis title of the boxplot (e.g.
#'   "Interaction", used for the axis text specifically).
#'
#' @return A named list with two ggplot objects: `scatter_plot` (males vs
#'   females mean log2FC per gene) and `box_plot` (per-sample projection
#'   coefficient by group/Condition).
sexLogFC_stressSignature_plts <- function(SE, gene_list, text1, text2){

  stressed_animals <- list(
  males = SE$SampleName[SE$Sex == "Male" & SE$Condition == "Stress"],
  females = SE$SampleName[SE$Sex == "Female" & SE$Condition == "Stress"])

  # 1. scatter plot
  log2fc_mat <- assay(SE, "log2FC")  # gene x sample matrix
  rownames(log2fc_mat) <- sub(".*\\.", "", rownames(log2fc_mat))
  gene_list_renamed <- sub(".*\\.", "", gene_list)
  
  x_vals <- log2fc_mat[gene_list_renamed, stressed_animals$males]
  y_vals <- log2fc_mat[gene_list_renamed, stressed_animals$females]
  
  x_means <- rowMeans(x_vals, na.rm = TRUE)
  y_means <- rowMeans(y_vals, na.rm = TRUE)
  
  df <- data.frame(males = x_means, females = y_means)
    
  # A label is added if the point is among the most significant to avoid cluttering and label interesting genes

  df$label <- ifelse(rownames(df) %in% head(gene_list_renamed, 20) | 
                       (abs(df$males) > 0.7 | abs(df$females) > 0.7), rownames(df), "")
  
  p1 <- ggplot(df, aes(x = males, y = females)) +
    geom_point(size = 0.7) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "gray40") +
    geom_hline(yintercept = 0, color = "gray40") +
    xlab("log2FC (Males)") +
    ylab("log2FC (Females)") +
    ggtitle(sprintf("%s DEGs: Males vs Females", text1)) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = Inf, force_pull = 5, min.segment.length = 0) +
    theme_bw() +
    geom_smooth(method = MASS::rlm)
  
  # 2. boxplot
  log2fc_mat <- assay(SE, "log2FC")
  e <- log2fc_mat[gene_list, ] # subset the SV corrected log2FC matrix to the stagewise corrected stress DEGs
  stressed_samples <- unlist(stressed_animals)
  x <- rowMeans(e[, stressed_samples])
  SE$coef2 <- sapply(colnames(e), function(current_sample) {
    y <- e[, current_sample]
    coef(summary(MASS::rlm(y ~ 0 + x)))[1, 1]
  })
  
  d <- as.data.frame(colData(SE))
  
  p2 <- ggplot(d, aes(group, coef2, fill = Condition)) +
    geom_boxplot() +
    theme_bw() +    
    ggtitle(sprintf("%s DEGs signature plot", text1)) +
    geom_point(position = position_jitterdodge(), size = 1) +
    ylab(sprintf("%s signature projection coefficient", text2)) + 
    stat_compare_means(
    comparisons = list(
      c("Male_Ctrl", "Female_Ctrl"),
      c("Male_Stress", "Female_Stress")
    ),
    method = "t.test")
  
  return(list(scatter_plot = p1, box_plot = p2))
}

#' Plot a volcano plot for the differential analysis results
#'
#' @param df the result table from differential analysis
#' @param logFC logFC value (from the result table in our case)
#' @param signif FDR value (from the result table in our case)
#' @param labels label special genes
#' @param sig.th threshold of the FDR
#' @param lfc.th threshold of logFC
#'
#' @return a volcano plot
volcano <- function(df, logFC, signif, labels=FALSE, sig.th=0.05, lfc.th=0.1){
  df <- dplyr::mutate(df, type=case_when(
        {{ signif }} < sig.th & {{ logFC }} > lfc.th ~ "Upregulated",
        {{ signif }} < sig.th & {{ logFC }} < -lfc.th ~ "Downregulated",
        TRUE ~ "n.s."
      ))
  if(isTRUE(labels)){
    df$label <- gsub(".+\\.","",row.names(df))
    df$label[which(df$type=="n.s.")] <- NA_character_
    labels <- "label"
  }
  if(is.character(labels) && (length(labels)>1 ||
      all(labels %in% row.names(df)))){
    df$label <- gsub(".+\\.","",row.names(df))
    df[!(row.names(df) %in% labels), "label"] <- NA_character_
    labels <- "label"
  }
  
  p <- ggplot(df, aes(
    x = {{logFC}},
    y = -log10({{signif}}),
    color = type)
  ) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "n.s." = "grey10"),
      name = ""
    ) + 
    geom_vline(xintercept = -lfc.th, linetype = "dashed", colour = "darkgray") +
    geom_vline(xintercept = lfc.th, linetype = "dashed", colour = "darkgray") +
    geom_hline(yintercept = -log10(sig.th), linetype = "dashed", colour = "darkgray") +
    theme_bw()

  if(!isFALSE(labels)){
    p <- p + geom_text_repel(aes(label = label), size = 2, max.overlaps = Inf, 
                             box.padding = 0.2, point.padding = 0.1, show.legend = FALSE)
  }

  return(p)
}

#' Two-stage FDR correction wrapper (stageR)
#'
#' Applies the stageR two-stage testing procedure to control the FDR when
#' a "screening" test is used to select genes for a more specific "confirmation" test. 
#' Testing the confirmation p-value only among screen-significant genes, 
#' with a naive FDR correction on the confirmation p-values alone, would 
#' inflate the false discovery rate; stageR corrects for this two-step selection process.
#'
#' @param screen A data.frame of screening results (e.g. from
#'   `topTags()`), with row names = gene IDs and a `FDR` column already
#'   containing (screen-level) adjusted p-values.
#' @param conf A data.frame of confirmation results (e.g. the interaction
#'   contrast), with row names = gene IDs and a `PValue` column of raw
#'   (unadjusted) p-values. 
#' @param alpha Nominal significance level passed to `stageWiseAdjustment()`.
#'
#' @return `conf`, with an added `staged.FDR` column (the stage-wise
#'   adjusted p-value for the confirmation test), sorted by `staged.FDR`
#'   ascending.
#' 
stagedFDRwrapper <- function(screen, conf, alpha=0.05){
  staged <- stageR(pScreen=screen$FDR, pConfirmation=cbind(interaction=conf[row.names(screen),"PValue"]), pScreenAdjusted = TRUE)
  staged <- stageWiseAdjustment(staged, method="none", alpha=alpha)
  stage.padj <- getAdjustedPValues(staged, onlySignificantGenes=FALSE, order=FALSE)
  row.names(stage.padj) <- row.names(screen)
  conf$staged.FDR <- as.data.frame(stage.padj)[row.names(conf), "interaction"]
  conf[order(conf$staged.FDR),]
}


#' Approximate logFC standard error and confidence bounds from DEA result
addlogFCSE <- function(x){
  x$SE <- abs(x$logFC) / sqrt(x$`F`)
  x$high <- x$logFC+x$SE
  x$low <- x$logFC-x$SE
  x
}

#' Extract log2FC data and create a dataframe for plotting
get_logFC_df <- function(SE, gene_list, dataset_name) {
  # Get stressed animal sample names
  stressed_animals <- list(
    males = SE$SampleName[SE$Sex == "Male" & SE$Condition == "Stress"],
    females = SE$SampleName[SE$Sex == "Female" & SE$Condition == "Stress"]
  )
  
  # Get log2FC matrix and clean up rownames
  log2fc_mat <- assay(SE, "log2FC")
  rownames(log2fc_mat) <- clean_name(rownames(log2fc_mat))
  gene_list_renamed <- clean_name(gene_list)
  
  # Ensure that the genes exist in the matrix
  common_genes <- intersect(gene_list_renamed, rownames(log2fc_mat))
  
  # Get log2FC values
  x_vals <- log2fc_mat[common_genes, stressed_animals$males]
  y_vals <- log2fc_mat[common_genes, stressed_animals$females]
  
  # Calculate row means
  x_means <- rowMeans(x_vals, na.rm = TRUE)
  y_means <- rowMeans(y_vals, na.rm = TRUE)
  
  # Create a dataframe
  df <- data.frame(
    males = x_means,
    females = y_means,
    dataset = dataset_name,
    gene = common_genes
  )
  
  return(df)
}


#' plot the log2FC of certain genes
genePlot <- function(SE, genes, assay = "log2FC", ...){
  d <- meltSE(SE, genes)
  p <- ggplot(d, aes(Condition, .data[[assay]], color=Sex)) + 
          stat_summary(geom="crossbar", position = position_dodge2(width=0.8)) + 
          geom_beeswarm(dodge.width = 0.8) + theme_classic() +
          theme(strip.background = element_blank(), strip.text = element_text(size=10))
 if(length(unique(d$feature))>1){
   p <- p + facet_wrap(~feature, ...)
 }else{
   p + ylab(paste0(unique(d$feature), " ", assay))
 }
  p
}

#' remove the NA rows and columns from RRHO result matrix, to plot better heatmap. And rotate the matrix to the correct form
make_better_RRHO_matrix <- function(mat){
  na_rows <- which(apply(mat, 1, function(x) all(is.na(x))))
  na_cols <- which(apply(mat, 2, function(x) all(is.na(x))))
  # keep only 1 NA row/col as a thin gap instead of all of them
  if(length(na_rows) > 1) mat <- mat[-na_rows[-1], ]
  if(length(na_cols) > 1) mat <- mat[, -na_cols[-1]]
  mat <- t(mat[,ncol(mat):1])
  return(mat)
}

#' plot RRHO heatmap manually
RRHO_plot <- function(mat, colorFun, col_title, row_title, name =  "-log10(P-value)", fontsize = 10, show_legend = TRUE, ...){
  Heatmap(mat, 
          cluster_rows = FALSE, 
          cluster_columns = FALSE, 
          na_col = "white", 
          col = colorFun,
          name = name,
          use_raster = TRUE,
          show_heatmap_legend = show_legend,
          column_title = col_title,
          column_title_side = "bottom",
          row_title = row_title,
          column_title_gp = gpar(fontsize = fontsize),
          row_title_gp = gpar(fontsize = fontsize),
          ...)
}

#' Extract the legends of a heatmap
get_ht_legend <- function(ht, gap = unit(4, "mm")) {
   # annotation legends (Condition, Sex, ...)
   ann <- lapply(ht@top_annotation@anno_list, function(a) {
     if (is.null(a@color_mapping)) return(NULL)
     do.call(color_mapping_legend,
             c(list(a@color_mapping, plot = FALSE), a@legend_param))
   })
   # heatmap (z-score) legend, with your heatmap_legend_param (title, at, ...)
   hm <- do.call(color_mapping_legend,
                 c(list(ht@matrix_color_mapping, plot = FALSE), ht@matrix_legend_param))
 
   do.call(packLegend, c(Filter(Negate(is.null), ann), list(hm),
                         list(direction = "vertical", gap = gap)))
 }

 
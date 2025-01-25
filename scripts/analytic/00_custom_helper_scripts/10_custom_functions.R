
pipe.fast.seu <- function(seu) {
  seu <- ScaleData(seu, do.scale = F, do.center = F)
  seu <- RunPCA(seu, features = rownames(seu))
  seu <- FindNeighbors(object = seu)
  seu <- FindClusters(object = seu)
  seu <- RunUMAP(seu, dims = 1:30)
  return(seu)
}

#Custom functions to convert between beta and mvalues
m_value_to_beta_value <- function(mvalue) {
  beta <- 2^mvalue/(2^mvalue + 1)
  return(beta)
}

beta_value_to_m_value <- function(beta) {
  mvalue <- log2(beta/(1-beta))
  return(mvalue)
}

#Here we can define another function to use to look at principle component plots by different attributes.
pcplot <- function(pheno){
  ggpairs(dat[, c("Dim.1", "Dim.2", "Dim.3", "Dim.4", "Dim.5")], aes(color = factor(dat %>% pull(pheno))))
}

# covariates.factor are converted to factors
pc_covar_test <- 
  function(beta, pd, covariates, covariates.factor) {
    set.seed(982278)
    pca.dat = beta
    dim(pca.dat)
    pca.dat = na.omit(pca.dat)
    dim(pca.dat)
    pca.dat = t(pca.dat)
    
    res.pca <- PCA(pca.dat, scale.unit = FALSE, ncp = 15, graph = FALSE)
    
    print(fviz_eig(res.pca, addlabels = TRUE))#, ylim = c(0, 100))
    
    var <- get_pca_var(res.pca)
    
    #pc.coord <- var$coord %>% as.data.frame %>% rownames_to_column(var = "code")
    # Pull individual coordinates
    pc.coord <- res.pca$ind$coord %>% as.data.frame
    pc.coord$Sample_Name <- pd$Sample_Name
    
    dat <- left_join(pd, pc.coord, by = "Sample_Name")
    
    pcs.keep.i <- grepl("Dim", names(dat))
    pcs.keep <- names(dat[pcs.keep.i])
    
    dat <- 
      dat %>%
      mutate(across(.cols = matches(covariates.factor), .fns = factor))
    
    map_outcome_to_predictors <- 
      function(outcome, predictors) {
        n.predictors = length(predictors)
        outcome.vector = rep(outcome, times = n.predictors)
        strings = paste0(outcome.vector, " ~ ", predictors)
        formulas <- map(strings, as.formula)
        return(formulas)
      }
    
    formulas = map_outcome_to_predictors(pcs.keep, covariates)
    
    models = map(formulas, ~ lm(.x, data = dat))
    
    models.stats <- map(models, glance)
    
    names(models.stats) <- as.character(formulas)
    
    model.dat <- lapply(models.stats, unlist) %>% as.data.frame() %>% t() %>% as.data.frame
    model.dat <- rownames_to_column(model.dat, var = "formula")
    
    outcome <- str_split(model.dat$formula, pattern = "\\.\\.\\.", simplify = TRUE)[,1]
    predictor <- str_split(model.dat$formula, pattern = "\\.\\.\\.", simplify = TRUE)[,2]
    
    model.dat$predictor <- predictor
    model.dat$outcome <- sub("Dim\\.", "PC", outcome)
    
    model.dat <- 
      model.dat %>%
      mutate(
        pval.cutoff =
          case_when(
            p.value.value > 0.1 ~ "Insignificant",
            p.value.value < 0.001 ~ "p < 0.001",
            p.value.value < 0.01 ~ "p < 0.01",
            p.value.value < 0.1 ~ "p < 0.1"
          )
      )
    
    model.dat$pval.cutoff <- fct_relevel(model.dat$pval.cutoff, "Insignificant", "p < 0.1", "p < 0.01", "p < 0.001")
    
    ggplot(data = model.dat %>% filter(outcome %in% c("PC1", "PC2", "PC3", "PC4", "PC5")),
           mapping = aes(
             x = outcome,
             y = predictor,
             fill = pval.cutoff
           )) +
      geom_tile() +
      xlab("") +
      ylab("") +
      theme_bw() +
      scale_fill_grey(start = 0.9, end = 0.1)
  }

run_limma_model <- function(expr, pd) {
  
  # Construct the expression set object required to run limma
  # here can use mvalue or beta value; we can use mvalue for statistical testing but beta value cutoffs for presentation
  # Standard approach is to test and present on beta value scale
  d <- ExpressionSet(expr)
  pData(d) <- pd
  
  #Fit the linear model
  #design <- model.matrix(~ 0 + phenotype, data = d)
  design <- model.matrix(~ 0 + phenotype + Placenta, data = d)
  colnames(design) <- make.names(colnames(design))
  colnames(design)
  fit <- lmFit(d, design)
  fit <- eBayes(fit, robust = TRUE)
  
  #Identify the contrasts you would like to investigate
  contrast.matrix <- makeContrasts(phenotypeCytotrophoblast - (phenotypeEndothelial + phenotypeHofbauer + phenotypeStromal + phenotypeSyncytiotrophoblast)/4,
                                   phenotypeEndothelial - (phenotypeCytotrophoblast + phenotypeHofbauer + phenotypeStromal + phenotypeSyncytiotrophoblast)/4,
                                   phenotypeHofbauer - (phenotypeCytotrophoblast + phenotypeEndothelial + phenotypeStromal + phenotypeSyncytiotrophoblast)/4,
                                   phenotypeStromal - (phenotypeCytotrophoblast + phenotypeEndothelial + phenotypeHofbauer + phenotypeSyncytiotrophoblast)/4,
                                   phenotypeSyncytiotrophoblast - (phenotypeCytotrophoblast + phenotypeEndothelial + phenotypeHofbauer + phenotypeStromal)/4,
                                   levels = design)
  colnames(contrast.matrix) <- c("Cytotrophoblast", "Endothelial", "Hofbauer", "Stromal", "Syncytiotrophoblast")
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2, robust = TRUE)
  
  return(fit2)
  
}

run_limma_model_450k_all_cell_types <- function(expr, pd) {
  
  # Construct the expression set object required to run limma
  # here can use mvalue or beta value; we can use mvalue for statistical testing but beta value cutoffs for presentation
  # Standard approach is to test and present on beta value scale
  d <- ExpressionSet(expr)
  pData(d) <- pd
  
  #Fit the linear model
  #design <- model.matrix(~ 0 + phenotype, data = d)
  
  # Blocking by donor not possible for the immune and nRBC cell types
  design <- model.matrix(~ 0 + CellType + Sex, data = d)
  colnames(design) <- make.names(colnames(design))
  colnames(design)
  fit <- lmFit(d, design)
  fit <- eBayes(fit, robust = TRUE)
  
  colnames(d)
  
  # Get CellType groups of interest
  groups <- colnames(design)[grep("CellType", colnames(design))]
  #groups.placenta
  
  # Get the names of contrasts
  contrasts <- map(
    groups,
    ~ paste0(.x, " - (", paste(groups[groups != .x], collapse = " + "), ")/", length(groups[groups != .x]))
  ) %>%
    unlist()
  
  contrast.matrix <- makeContrasts(Cytotrophoblast = CellTypeCytotrophoblast - (CellTypeSyncytiotrophoblast + CellTypeHofbauer + CellTypeStromal + CellTypeEndothelial + CellTypeNeu + CellTypeNK + CellTypeBcell + CellTypeCD4T + CellTypeCD8T + CellTypeMono + CellTypenRBC)/11,
                                   Syncytiotrophoblast = CellTypeSyncytiotrophoblast - (CellTypeCytotrophoblast + CellTypeHofbauer + CellTypeStromal + CellTypeEndothelial + CellTypeNeu + CellTypeNK + CellTypeBcell + CellTypeCD4T + CellTypeCD8T + CellTypeMono + CellTypenRBC)/11,
                                   Hofbauer = CellTypeHofbauer - (CellTypeCytotrophoblast + CellTypeSyncytiotrophoblast + CellTypeStromal + CellTypeEndothelial + CellTypeNeu + CellTypeNK + CellTypeBcell + CellTypeCD4T + CellTypeCD8T + CellTypeMono + CellTypenRBC)/11,
                                   Stromal = CellTypeStromal - (CellTypeCytotrophoblast + CellTypeSyncytiotrophoblast + CellTypeHofbauer + CellTypeEndothelial + CellTypeNeu + CellTypeNK + CellTypeBcell + CellTypeCD4T + CellTypeCD8T + CellTypeMono + CellTypenRBC)/11,
                                   Endothelial = CellTypeEndothelial - (CellTypeCytotrophoblast + CellTypeSyncytiotrophoblast + CellTypeHofbauer + CellTypeStromal + CellTypeNeu + CellTypeNK + CellTypeBcell + CellTypeCD4T + CellTypeCD8T + CellTypeMono + CellTypenRBC)/11,
                                   Neutrophil = CellTypeNeu - (CellTypeCytotrophoblast + CellTypeSyncytiotrophoblast + CellTypeHofbauer + CellTypeStromal + CellTypeEndothelial + CellTypeNK + CellTypeBcell + CellTypeCD4T + CellTypeCD8T + CellTypeMono + CellTypenRBC)/11,
                                   NaturalKiller = CellTypeNK - (CellTypeCytotrophoblast + CellTypeSyncytiotrophoblast + CellTypeHofbauer + CellTypeStromal + CellTypeEndothelial + CellTypeNeu + CellTypeBcell + CellTypeCD4T + CellTypeCD8T + CellTypeMono + CellTypenRBC)/11,
                                   BCell = CellTypeBcell - (CellTypeCytotrophoblast + CellTypeSyncytiotrophoblast + CellTypeHofbauer + CellTypeStromal + CellTypeEndothelial + CellTypeNeu + CellTypeNK + CellTypeCD4T + CellTypeCD8T + CellTypeMono + CellTypenRBC)/11,
                                   CD4T = CellTypeCD4T - (CellTypeCytotrophoblast + CellTypeSyncytiotrophoblast + CellTypeHofbauer + CellTypeStromal + CellTypeEndothelial + CellTypeNeu + CellTypeNK + CellTypeBcell + CellTypeCD8T + CellTypeMono + CellTypenRBC)/11,
                                   CD8T = CellTypeCD8T - (CellTypeCytotrophoblast + CellTypeSyncytiotrophoblast + CellTypeHofbauer + CellTypeStromal + CellTypeEndothelial + CellTypeNeu + CellTypeNK + CellTypeBcell + CellTypeCD4T + CellTypeMono + CellTypenRBC)/11,
                                   Monocyte = CellTypeMono - (CellTypeCytotrophoblast + CellTypeSyncytiotrophoblast + CellTypeHofbauer + CellTypeStromal + CellTypeEndothelial + CellTypeNeu + CellTypeNK + CellTypeBcell + CellTypeCD4T + CellTypeCD8T + CellTypenRBC)/11,
                                   nRBC = CellTypenRBC - (CellTypeCytotrophoblast + CellTypeSyncytiotrophoblast + CellTypeHofbauer + CellTypeStromal + CellTypeEndothelial + CellTypeNeu + CellTypeNK + CellTypeBcell + CellTypeCD4T + CellTypeCD8T + CellTypeMono)/11,
                                   levels = design)

  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2, robust = TRUE)
  
  return(fit2)
  
}

graph_limma_model <-
  function(contrast.name,
           hits,
           FDR_threshold,
           logFC_threshold,
           genes.of.interest,
           ...) {
    #Histogram of raw p-values
    hist(
      hits$P.Value,
      col = brewer.pal(3, name = "Set2")[1],
      main = contrast.name,
      xlab = "P-values"
    )
    
    #Histogram of FDR-adjusted p-values.
    hist(
      hits$adj.P.Val,
      col = brewer.pal(3, name = "Set2")[1],
      main = contrast.name,
      xlab = "FDR-adjusted p-values"
    )
    
    # Print summary of hits
    hits$sig <- (hits$adj.P.Val < FDR_threshold) & (abs(hits$logFC) > logFC_threshold)
    print(paste0(
      length(which(hits$sig)),
      " significant hits at FDR-adjusted cutoff of ",
      FDR_threshold
    ))
    
    # Label factor levels for legend purposes
    hits$sig <- factor(hits$sig,
                       levels = c(FALSE, TRUE),
                       labels = c("Does not meet criteria", paste0("Meets FDR-adjusted p-value < ", FDR_threshold, " and logFC > ", logFC_threshold)))
    
    # Move gene symbols to column
    res <- rownames_to_column(hits, "gene_symbol")
    # Get the top hypomethylated hits
    downreg <- res %>%
      filter(adj.P.Val < FDR_threshold) %>%
      filter(logFC < 0) %>%
      filter(abs(logFC) > logFC_threshold) %>%
      arrange(logFC)
    #Get the top hypermethylated hits.
    upreg <- res %>%
      filter(adj.P.Val < FDR_threshold) %>%
      filter(logFC > 0) %>%
      filter(abs(logFC) > logFC_threshold) %>%
      arrange(desc(logFC))
    #Label genes of interest by name, the top 10 hypermethylated and hypometylated hits and hand-chosen genes of interest.
    num.to.plot <- 10
    top.upreg <- upreg$gene_symbol[1:num.to.plot]
    top.downreg <- downreg$gene_symbol[1:num.to.plot]
    genes.to.plot <- c(top.upreg, top.downreg, genes.of.interest)
    res$label[(res$gene_symbol %in% genes.to.plot)] <- T
    # Volcano plot
    plot <- ggplot(res) +
      geom_point(aes(
        x = logFC,
        y = -log10(adj.P.Val),
        colour = sig
      )) +
      geom_text_repel(aes(
        x = logFC,
        y = -log10(adj.P.Val),
        label = ifelse(label, gene_symbol, "")
      )) +
      geom_hline(yintercept = -log10(FDR_threshold),
                 linetype = "dotted") +
      geom_vline(xintercept = -logFC_threshold,
                 linetype = "dotted") +
      geom_vline(xintercept = logFC_threshold,
                 linetype = "dotted") +
      ggtitle(contrast.name) +
      xlab("logFC of  M-value") +
      ylab("-log10 adjusted p-value") +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))
      )
    
    # Create a label dataframe for geom_label annotation
    label.df <- data.frame(
      label = paste0(
        dim(upreg)[1],
        " hypermethylated sites\n",
        paste0(dim(downreg)[1], " hypomethylated sites")
      ),
      x = Inf,
      y = Inf,
      hjustvar = 1,
      vjustvar = 1
    )
    
    # Add geom_label annotation
    plot <-
      plot + geom_label(
        data = label.df,
        aes(
          x = x,
          y = y,
          label = label,
          hjust = hjustvar,
          vjust = vjustvar
        ),
        #size = 6,
        label.padding = unit(0.25, "lines")
      )
    
    # Print plot to output
    print(plot)
    # Create list of volcano plot and models result for export
    return(list
           (volcano = plot,
             res = res))
  }


graph_limma_model_scattermore_beta <-
  function(contrast.name,
           hits,
           FDR_threshold,
           beta_threshold,
           genes.of.interest,
           ...) {
    #Histogram of raw p-values
    hist(
      hits$P.Value,
      col = brewer.pal(3, name = "Set2")[1],
      main = contrast.name,
      xlab = "P-values"
    )
    
    #Histogram of FDR-adjusted p-values.
    hist(
      hits$adj.P.Val,
      col = brewer.pal(3, name = "Set2")[1],
      main = contrast.name,
      xlab = "FDR-adjusted p-values"
    )
    
    # Print summary of hits
    hits$sig <- (hits$adj.P.Val < FDR_threshold) & (abs(hits$logFC) > beta_threshold)
    print(paste0(
      length(which(hits$sig)),
      " significant hits at FDR-adjusted cutoff of ",
      FDR_threshold
      , " and log2-fold threshold of ", beta_threshold, "."))
    
    # Label factor levels for legend purposes
    hits$sig <- factor(hits$sig,
                       levels = c(FALSE, TRUE),
                       labels = c("Does not meet criteria", paste0("Meets FDR-adjusted p-value < ", FDR_threshold, " and methylation rate difference > ", beta_threshold*100, "%")))
    sig.label <- paste0("Meets FDR-adjusted p-value < ", FDR_threshold, " and methylation rate difference > ", beta_threshold*100, "%")
    
    # Move gene symbols to column
    res <- rownames_to_column(hits, "gene_symbol")
    # Get the top hypomethylated hits
    downreg <- res %>%
      filter(adj.P.Val < FDR_threshold) %>%
      filter(logFC < 0) %>%
      filter(abs(logFC) > beta_threshold) %>%
      arrange(logFC)
    #Get the top hypermethylated hits.
    upreg <- res %>%
      filter(adj.P.Val < FDR_threshold) %>%
      filter(logFC > 0) %>%
      filter(abs(logFC) > beta_threshold) %>%
      arrange(desc(logFC))
    #Label genes of interest by name, the top 10 hypermethylated and hypometylated hits and hand-chosen genes of interest.
    num.to.plot <- 10
    top.upreg <- upreg$gene_symbol[1:num.to.plot]
    top.downreg <- downreg$gene_symbol[1:num.to.plot]
    genes.to.plot <- c(top.upreg, top.downreg, genes.of.interest)
    res$label[(res$gene_symbol %in% genes.to.plot)] <- T
    # Volcano plot
    plot <- ggplot(res) +
      geom_scattermore(aes(
        x = logFC,
        y = -log10(adj.P.Val),
        colour = sig
      )) +
      geom_text_repel(aes(
        x = logFC,
        y = -log10(adj.P.Val),
        label = ifelse(label, gene_symbol, "")
      )) +
      geom_hline(yintercept = -log10(FDR_threshold),
                 linetype = "dotted") +
      geom_vline(xintercept = -beta_threshold,
                 linetype = "dotted") +
      geom_vline(xintercept = beta_threshold,
                 linetype = "dotted") +
      theme_bw() +
      # Manually set colors
      scale_color_manual(values = c("gray", "red")) +
      ggtitle(paste0(contrast.name, " vs. others")) +
      scale_x_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 10)) +
      xlab("DNA Methylation Difference") +
      ylab("-log10 adjusted p-value") +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))
      )
    
    # Create a label dataframe for geom_label annotation
    label.df <- data.frame(
      label = paste0(
        dim(upreg)[1],
        " hypermethylated sites\n",
        paste0(dim(downreg)[1], " hypomethylated sites")
      ),
      x = -Inf,
      y = -Inf,
      hjustvar = 0,
      vjustvar = 0
    )
    
    # Add geom_label annotation
    plot <-
      plot + geom_label(
        data = label.df,
        aes(
          x = x,
          y = y,
          label = label,
          hjust = hjustvar,
          vjust = vjustvar
        ),
        #size = 6,
        label.padding = unit(0.25, "lines")
      )
    
    # Print plot to output
    print(plot)
    # Create list of volcano plot and models result for export
    return(list
           (volcano = plot,
             res = res))
  }

graph_limma_model_scattermore <-
  function(contrast.name,
           hits,
           FDR_threshold,
           logFC_threshold,
           genes.of.interest,
           ...) {
    #Histogram of raw p-values
    hist(
      hits$P.Value,
      col = brewer.pal(3, name = "Set2")[1],
      main = contrast.name,
      xlab = "P-values"
    )
    
    #Histogram of FDR-adjusted p-values.
    hist(
      hits$adj.P.Val,
      col = brewer.pal(3, name = "Set2")[1],
      main = contrast.name,
      xlab = "FDR-adjusted p-values"
    )
    
    # Print summary of hits
    hits$sig <- (hits$adj.P.Val < FDR_threshold) & (abs(hits$logFC) > logFC_threshold)
    print(paste0(
      length(which(hits$sig)),
      " significant hits at FDR-adjusted cutoff of ",
      FDR_threshold
    , " and log2-fold threshold of ", logFC_threshold, "."))
    
    # Label factor levels for legend purposes
    hits$sig <- factor(hits$sig,
                       levels = c(FALSE, TRUE),
                       labels = c("Does not meet criteria", paste0("Meets FDR-adjusted p-value < ", FDR_threshold, " and log2FC > ", logFC_threshold)))
    sig.label <- paste0("Meets FDR-adjusted p-value < ", FDR_threshold, " and log2FC > ", logFC_threshold)
    
    # Move gene symbols to column
    res <- rownames_to_column(hits, "gene_symbol")
    # Get the top hypomethylated hits
    downreg <- res %>%
      filter(adj.P.Val < FDR_threshold) %>%
      filter(logFC < 0) %>%
      filter(abs(logFC) > logFC_threshold) %>%
      arrange(logFC)
    #Get the top hypermethylated hits.
    upreg <- res %>%
      filter(adj.P.Val < FDR_threshold) %>%
      filter(logFC > 0) %>%
      filter(abs(logFC) > logFC_threshold) %>%
      arrange(desc(logFC))
    #Label genes of interest by name, the top 10 hypermethylated and hypometylated hits and hand-chosen genes of interest.
    num.to.plot <- 10
    top.upreg <- upreg$gene_symbol[1:num.to.plot]
    top.downreg <- downreg$gene_symbol[1:num.to.plot]
    genes.to.plot <- c(top.upreg, top.downreg, genes.of.interest)
    res$label[(res$gene_symbol %in% genes.to.plot)] <- T
    # Volcano plot
    plot <- ggplot(res) +
      geom_scattermore(aes(
        x = logFC,
        y = -log10(adj.P.Val),
        colour = sig
      )) +
      geom_text_repel(aes(
        x = logFC,
        y = -log10(adj.P.Val),
        label = ifelse(label, gene_symbol, "")
      )) +
      geom_hline(yintercept = -log10(FDR_threshold),
                 linetype = "dotted") +
      geom_vline(xintercept = -logFC_threshold,
                 linetype = "dotted") +
      geom_vline(xintercept = logFC_threshold,
                 linetype = "dotted") +
      # Manually set colors
      scale_color_manual(
        values = c("Does not meet criteria" = "gray",
                   `sig.label` = "salmon")) +
      ggtitle(paste0(contrast.name, " vs. others")) +
      xlab("log2-fold change of  M-value") +
      ylab("-log10 adjusted p-value") +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))
      )
    
    # Create a label dataframe for geom_label annotation
    label.df <- data.frame(
      label = paste0(
        dim(upreg)[1],
        " hypermethylated sites\n",
        paste0(dim(downreg)[1], " hypomethylated sites")
      ),
      x = -Inf,
      y = -Inf,
      hjustvar = 0,
      vjustvar = 0
    )
    
    # Add geom_label annotation
    plot <-
      plot + geom_label(
        data = label.df,
        aes(
          x = x,
          y = y,
          label = label,
          hjust = hjustvar,
          vjust = vjustvar
        ),
        #size = 6,
        label.padding = unit(0.25, "lines")
      )
    
    # Print plot to output
    print(plot)
    # Create list of volcano plot and models result for export
    return(list
           (volcano = plot,
             res = res))
  }

# Taken from Hopkins cluster DNAm code
graph_limma_model_p_val_raw <-
  function(contrast.name,
           hits,
           P.Value_threshold,
           logFC_threshold,
           genes.of.interest,
           ...) {
    #Histogram of raw p-values
    hist(
      hits$P.Value,
      col = brewer.pal(3, name = "Set2")[1],
      main = contrast.name,
      xlab = "P-values"
    )
    
    #Histogram of p-value p-values.
    hist(
      hits$P.Value,
      col = brewer.pal(3, name = "Set2")[1],
      main = contrast.name,
      xlab = "p-value p-values"
    )
    
    # Print summary of hits
    hits$sig <- (hits$P.Value < P.Value_threshold) & (abs(hits$logFC) > logFC_threshold)
    print(paste0(
      length(which(hits$sig)),
      " significant hits at p-value cutoff of ",
      P.Value_threshold
    ))
    
    # Label factor levels for legend purposes
    hits$sig <- factor(hits$sig,
                       levels = c(FALSE, TRUE),
                       labels = c("Does not meet criteria", paste0("Meets p-value p-value < ", P.Value_threshold, " and absolute percentage difference > ", logFC_threshold)))
    
    # Move gene symbols to column
    res <- rownames_to_column(hits, "gene_symbol")
    # Get the top hypomethylated hits
    downreg <- res %>%
      filter(P.Value < P.Value_threshold) %>%
      filter(logFC < 0) %>%
      filter(abs(logFC) > logFC_threshold) %>%
      arrange(logFC)
    #Get the top hypermethylated hits.
    upreg <- res %>%
      filter(P.Value < P.Value_threshold) %>%
      filter(logFC > 0) %>%
      filter(abs(logFC) > logFC_threshold) %>%
      arrange(desc(logFC))
    #Label genes of interest by name, the top 10 hypermethylated and hypometylated hits and hand-chosen sites of interest.
    num.to.plot <- 10
    top.upreg <- upreg$gene_symbol[1:num.to.plot]
    top.downreg <- downreg$gene_symbol[1:num.to.plot]
    genes.to.plot <- c(top.upreg, top.downreg, genes.of.interest)
    res$label[(res$gene_symbol %in% genes.to.plot)] <- T
    # Volcano plot
    plot <- ggplot(res) +
      geom_scattermore(aes(
        x = logFC,
        y = -log10(P.Value),
        colour = sig
      ),
      pixels = c(1080, 1080), pointsize = 3, alpha = 0.5) +
      #geom_text_repel(aes(
      #  x = logFC,
      #  y = -log10(P.Value),
      #  label = ifelse(label, gene_symbol, "")
      #)) +
      geom_hline(yintercept = -log10(1e-7),
                 linetype = "dotted") +
      geom_hline(yintercept = -log10(1e-2),
                linetype = "dotted") +
      geom_vline(xintercept = -logFC_threshold,
                 linetype = "dotted") +
      geom_vline(xintercept = logFC_threshold,
                 linetype = "dotted") +
      ggtitle(contrast.name) +
      xlab("Beta value difference") +
      ylab("-log10 p-value") +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))
      )
    
    # Create a label dataframe for geom_label annotation
    label.df <- data.frame(
      label = paste0(
        dim(upreg)[1],
        " hypermethylated sites\n",
        paste0(dim(downreg)[1], " hypomethylated sites")
      ),
      x = Inf,
      y = Inf,
      hjustvar = 1,
      vjustvar = 1
    )
    
    # Add geom_label annotation
    plot <-
      plot + geom_text(
        data = label.df,
        aes(
          x = x,
          y = y,
          label = label,
          hjust = hjustvar,
          vjust = vjustvar
        ),
        size = 12#,
        #label.padding = unit(0.25, "lines")
      )
    
    # Print plot to output
    #print(plot)
    # Create list of volcano plot and models result for export
    return(list
           (volcano = plot,
             res = res))
  }

#Custom function calculate M_model_M_mean delta beta-value method from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6449748/, for when covariate coefficients are unavailable (may be able to find them in the limma model)
M_model_M_mean <- function(m0, deltam) {
  delta_beta <- m_value_to_beta_value(m0 + deltam) - m_value_to_beta_value(m0)
  return(delta_beta)
}

#Custom function to pull hits from a limma fit and add beta.diff
get_limma_hits <- function(fit, coef) {
  
  # Pull m-value model results
  hits <- topTable(fit = fit, coef = coef, number = Inf, confint = TRUE)
  
  # Add beta value logFC to m-value model as beta.diff
  #hits$beta.diff <- M_model_M_mean(m0 = fit2$Amean, hits$logFC)
  
  # Subsample for graphing debugging purposes
  #hits<- hits[sample(nrow(hits), 0.01*nrow(hits)), ]
  
  return(hits)
}

#Custom function to pull hits from a limma fit and add beta.diff
get_limma_hits_beta <- function(fit, coef) {
  
  # Pull m-value model results
  hits <- topTable(fit = fit, coef = coef, number = Inf, confint = TRUE)

  # Subsample for graphing debugging purposes
  #hits<- hits[sample(nrow(hits), 0.01*nrow(hits)), ]
  
  return(hits)
}

#Version that graphs beta values
graph_limma_model_beta <-
  function(contrast.name,
           hits,
           FDR_threshold,
           beta.diff_threshold,
           genes.of.interest,
           ...) {
    
    # Convert add beta.diff as duplicate of LogFC column
    hits <- 
      hits %>%
      mutate(beta.diff = logFC)
    
    # Print summary of hits
    hits$sig <- (hits$adj.P.Val < FDR_threshold) & (abs(hits$beta.diff) > beta.diff_threshold)
    print(paste0(
      length(which(hits$sig)),
      " significant hits at p-value cutoff of ",
      FDR_threshold,
      " and absolute methylation difference > ",
      beta.diff_threshold*100,
      "%"
    ))
    
    # Label factor levels for legend purposes
    hits$sig <- factor(hits$sig,
                       levels = c(FALSE, TRUE),
                       labels = c("Does not meet criteria",
                                  paste0("Meets FDR-adjusted p-value < ",
                                         FDR_threshold,
                                         " and absolute methylation difference > ",
                                         beta.diff_threshold*100, "%")))
    
    # Move gene symbols to column
    res <- rownames_to_column(hits, "gene_symbol")
    
    # Get the top hypomethylated hits
    downreg <- res %>%
      filter(beta.diff < 0) %>%
      filter(sig == paste0("Meets FDR-adjusted p-value < ",
                           FDR_threshold,
                           " and absolute methylation difference > ",
                           beta.diff_threshold*100, "%")) %>%
      arrange(beta.diff)
    
    #Get the top hypermethylated hits.
    upreg <- res %>%
      filter(beta.diff > 0) %>%
      filter(sig == paste0("Meets FDR-adjusted p-value < ",
                           FDR_threshold,
                           " and absolute methylation difference > ",
                           beta.diff_threshold*100, "%")) %>%
      arrange(desc(beta.diff))
    
    #Label genes of interest by name, the top 10 hypermethylated and hypometylated hits and hand-chosen genes of interest.
    num.to.plot <- 10
    top.upreg <- upreg$gene_symbol[1:num.to.plot]
    top.downreg <- downreg$gene_symbol[1:num.to.plot]
    genes.to.plot <- c(top.upreg, top.downreg, genes.of.interest)
    res$label[(res$gene_symbol %in% genes.to.plot)] <- T
    
    # Volcano plot
    plot <- ggplot(res) +
      geom_scattermore(aes(
        x = beta.diff,
        y = -log10(adj.P.Val),
        colour = sig
      ), alpha = 0.25) +
      geom_text_repel(aes(
        x = beta.diff,
        y = -log10(adj.P.Val),
        label = ifelse(label, gene_symbol, "")
      ), alpha = 1) +
      geom_hline(yintercept = -log10(FDR_threshold),
                 linetype = "dotted") +
      geom_vline(xintercept = -beta.diff_threshold,
                 linetype = "dotted") +
      geom_vline(xintercept = beta.diff_threshold,
                 linetype = "dotted") +
      ggtitle(contrast.name) +
      xlab("Absolute average methylation difference") +
      ylab("-log10 adjusted p-value") +
      scale_x_continuous(labels = scales::percent) + 
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))
      )
    
    # Create a label dataframe for geom_label annotation
    label.df <- data.frame(
      label = paste0(
        dim(upreg)[1],
        " hypermethylated sites\n",
        paste0(dim(downreg)[1], " hypomethylated sites")
      ),
      x = -Inf,
      y = -Inf,
      hjustvar = 0,
      vjustvar = 0
    )
    
    # Add geom_label annotation
    plot <-
      plot + geom_label(data = label.df,
                        aes(
                          x = x,
                          y = y,
                          label = label,
                          hjust = hjustvar,
                          vjust = vjustvar
                        ),
                        label.padding = unit(0.25, "lines")
      )
    
    # Print plot to output
    #print(plot)
    # Create list of volcano plot and models result for export
    return(list
           (volcano = plot,
             res = res))
  }

graph_limma_model_m <-
  function(contrast.name,
           hits,
           FDR_threshold,
           beta.diff_threshold,
           genes.of.interest,
           ...) {
    
    # Print summary of hits
    hits$sig <- (hits$adj.P.Val < FDR_threshold) & (abs(hits$beta.diff) > beta.diff_threshold)
    print(paste0(
      length(which(hits$sig)),
      " significant hits at FDR-adjusted cutoff of ",
      FDR_threshold,
      " and beta difference threshold of ", beta.diff_threshold*100, "%"
    ))
    
    # Label factor levels for legend purposes
    hits$sig <- factor(hits$sig,
                       levels = c(FALSE, TRUE),
                       labels = c("Does not meet criteria", paste0("Meets FDR-adjusted p-value < ", FDR_threshold, " and site methylation difference > ", beta.diff_threshold*100, "%")))
    
    # Move gene symbols to column
    res <- rownames_to_column(hits, "gene_symbol")
    # Get the top hypomethylated hits
    downreg <- res %>%
      filter(adj.P.Val < FDR_threshold) %>%
      filter(beta.diff < 0) %>%
      filter(abs(beta.diff) > beta.diff_threshold) %>%
      arrange(beta.diff)
    #Get the top hypermethylated hits.
    upreg <- res %>%
      filter(adj.P.Val < FDR_threshold) %>%
      filter(beta.diff > 0) %>%
      filter(abs(beta.diff) > beta.diff_threshold) %>%
      arrange(desc(beta.diff))
    #Label genes of interest by name, the top 10 hypermethylated and hypometylated hits and hand-chosen genes of interest.
    num.to.plot <- 10
    top.upreg <- upreg$gene_symbol[1:num.to.plot]
    top.downreg <- downreg$gene_symbol[1:num.to.plot]
    genes.to.plot <- c(top.upreg, top.downreg, genes.of.interest)
    res$label[(res$gene_symbol %in% genes.to.plot)] <- T
    # Volcano plot
    plot <- ggplot(res) +
      geom_point(aes(
        x = logFC,
        y = -log10(adj.P.Val),
        colour = sig,
        shape = 1
      )) +
      scale_shape_identity() +
      geom_text_repel(aes(
        x = logFC,
        y = -log10(adj.P.Val),
        label = ifelse(label, gene_symbol, "")
      )) +
      geom_hline(yintercept = -log10(FDR_threshold),
                 linetype = "dotted") +
      #geom_vline(xintercept = -beta.diff_threshold,
      #           linetype = "dotted") +
      #geom_vline(xintercept = beta.diff_threshold,
      #           linetype = "dotted") +
      ggtitle(contrast.name) +
      xlab("log2 fold-change in methylation") +
      ylab("-log10 adjusted p-value") +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))
      )
    
    # Create a label dataframe for geom_label annotation
    label.df <- data.frame(
      label = paste0(
        dim(upreg)[1],
        " hypermethylated sites\n",
        paste0(dim(downreg)[1], " hypomethylated sites")
      ),
      x = -Inf,
      y = -Inf,
      hjustvar = 0,
      vjustvar = 0
    )
    
    # Add geom_label annotation
    plot <-
      plot + geom_label(data = label.df,
                        aes(
                          x = x,
                          y = y,
                          label = label,
                          hjust = hjustvar,
                          vjust = vjustvar
                        ),
                        label.padding = unit(0.25, "lines")
      )
    
    # Print plot to output
    print(plot)
    # Create list of volcano plot and models result for export
    return(list
           (volcano = plot,
             res = res))
  }

# Input ggplot to get x and y limits for graph
get_coord_limits <- function(ggplot) {
  
  # Get build options
  build <- ggplot_build(ggplot)
  # Get x-limits
  x <- build$layout$panel_scales_x[[1]]$range$range
  # Get y-limits
  y <- build$layout$panel_scales_y[[1]]$range$range
  
  # Return x and y limits as list
  return(list(
    "x" = x,
    "y" = y
  ))
}

# Get intersect_members for upsetR
get_intersect_members <- function (x, ...){
  require(dplyr)
  require(tibble)
  x <- x[,sapply(x, is.numeric)][,0<=colMeans(x[,sapply(x, is.numeric)],na.rm=T) & colMeans(x[,sapply(x, is.numeric)],na.rm=T)<=1]
  n <- names(x)
  x %>% rownames_to_column() -> x
  l <- c(...)
  a <- intersect(names(x), l)
  ar <- vector('list',length(n)+1)
  ar[[1]] <- x
  i=2
  for (item in n) {
    if (item %in% a){
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '>= 1')
        i <- i + 1
      }
    } else {
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '== 0')
        i <- i + 1
      }
    }
  }
  do.call(filter_, ar) %>% column_to_rownames() -> x
  return(x)
}

# Define a function to get top CpGs for input into missMethyl::gometh
get_cpg_top_hits_for_gometh <- function(res) {
  hits <-
    res %>%
    #dplyr::filter(logFC > 0) %>% # Limit to increased methylation difference; not applicable to methylation?
    filter(sig == levels(limma_results$Cytotrophoblast$sig)[2]) %>% # Limit to only CpGs that meet logFC and FDR cutoffs
    arrange(P.Value) %>% # Arrange by ascending p-value
    arrange(P.Value, desc(abs(logFC))) %>%
    slice_head(n = 10000) %>% # Pull top 10000 by p-value
    pull(gene_symbol) # Pull the cpg probe id for return
  return(hits)
}

# Define a function to get top CpGs for input into missMethyl::gometh
get_cpg_top_sig_hits_for_gometh_by_test_statistic <- function(res, top_n) {
  hits <-
    res %>%
    #dplyr::filter(logFC > 0) %>% # Limit to increased methylation difference; not applicable to methylation?
    filter(sig == levels(limma_results$Cytotrophoblast$sig)[2]) %>% # Limit to only CpGs that meet logFC and FDR cutoffs
    arrange(desc(abs(t))) %>% # Arrange by descending absolute value of test statistic
    slice_head(n = top_n) %>% # Pull top 10000 by p-value
    pull(gene_symbol) # Pull the cpg probe id for return
  return(hits)
}

# Define a function to get top cpgs per limma contrast
get_top_cpgs_per_limma_contrast <- function(res) {
  hits <-
    res %>%
    filter(sig == levels(limma_results$Cytotrophoblast$sig)[2]) %>% # Limit to only CpGs that meet logFC and FDR cutoffs
    arrange(desc(abs(`t`))) %>% # Arrange by moderated t-statistic
    slice_head(n = 100) %>% # Pull top 100 by moderated t-statistic
    pull(gene_symbol) # Pull the cpg probe id for return
  return(hits)
}

# Definate a function that accepts top cpgs as creates an ewastools reference matrix
generate_ewastools_ref_matrix <- function(top.cpgs) {
  
  beta.subset <- beta.celltypes[rownames(beta.celltypes) %in% top.cpgs,] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Sentrix.ID")
  
  beta.subset.anno <- 
    beta.subset %>%
    left_join(pd.celltypes %>%
                rownames_to_column(var = "Sentrix.ID") %>%
                dplyr::select(Sentrix.ID, phenotype)
    ) %>%
    dplyr::select(Sentrix.ID, phenotype, everything())
  
  #TODO there's a few missing beta values per top CpG
  beta.top.cpgs <-
    beta.subset.anno %>%
    pivot_longer(cols = all_of(top.cpgs), names_to = "Site", values_to = "Beta") %>%
    group_by(phenotype, Site) %>%
    summarise(Average = mean(Beta, na.rm = T)) %>%
    pivot_wider(names_from = phenotype, values_from = Average)
  
  ewastools.ref.matrix <- 
    beta.top.cpgs %>%
    column_to_rownames(var = "Site")
  
  return(ewastools.ref.matrix)
}


# Source code for validationCellType from minfi for local implementation
validationCellType <- function(Y, pheno, modelFix, modelBatch=NULL,
                               L.forFstat = NULL, verbose = FALSE){
  N <- dim(pheno)[1]
  pheno$y <- rep(0, N)
  xTest <- model.matrix(modelFix, pheno)
  sizeModel <- dim(xTest)[2]
  M <- dim(Y)[1]
  
  if (is.null(L.forFstat)) {
    # NOTE: All non-intercept coefficients
    L.forFstat <- diag(sizeModel)[-1,]
    colnames(L.forFstat) <- colnames(xTest)
    rownames(L.forFstat) <- colnames(xTest)[-1]
  }
  
  # Initialize various containers
  sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
  coefEsts <- matrix(NA, M, sizeModel)
  coefVcovs <- list()
  
  if (verbose) cat("[validationCellType] ")
  # Loop over each CpG
  for (j in seq_len(M)) {
    # Remove missing methylation values
    ii <- !is.na(Y[j, ])
    nObserved[j] <- sum(ii)
    pheno$y <- Y[j,]
    
    if (j %% round(M / 10) == 0 && verbose) cat(".") # Report progress
    
    # Try to fit a mixed model to adjust for plate
    try({
      if (!is.null(modelBatch)) {
        fit <- try(
          lme(modelFix, random = modelBatch, data = pheno[ii, ]))
        # NOTE: If LME can't be fit, just use OLS
        OLS <- inherits(fit, "try-error")
      } else {
        OLS <- TRUE
      }
      
      if (OLS) {
        fit <- lm(modelFix, data = pheno[ii, ])
        fitCoef <- fit$coef
        sigmaResid[j] <- summary(fit)$sigma
        sigmaIcept[j] <- 0
        nClusters[j] <- 0
      } else {
        fitCoef <- fit$coef$fixed
        sigmaResid[j] <- fit$sigma
        sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
        nClusters[j] <- length(fit$coef$random[[1]])
      }
      coefEsts[j,] <- fitCoef
      coefVcovs[[j]] <- vcov(fit)
      
      useCoef <- L.forFstat %*% fitCoef
      useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
      Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef)) / sizeModel
    })
  }
  if (verbose) cat(" done\n")
  
  # Name the rows so that they can be easily matched to the target data set
  rownames(coefEsts) <- rownames(Y)
  colnames(coefEsts) <- names(fitCoef)
  degFree <- nObserved - nClusters - sizeModel + 1
  
  # Get P values corresponding to F statistics
  Pval <- 1 - pf(Fstat, sizeModel, degFree)
  
  list(
    coefEsts = coefEsts,
    coefVcovs = coefVcovs,
    modelFix = modelFix,
    modelBatch = modelBatch,
    sigmaIcept = sigmaIcept,
    sigmaResid = sigmaResid,
    L.forFstat = L.forFstat,
    Pval = Pval,
    orderFstat = order(-Fstat),
    Fstat = Fstat,
    nClusters = nClusters,
    nObserved = nObserved,
    degFree = degFree)
}

# p is the input beta matrix
# pd.fn is the phenotype dataframe associated with that beta matrix
# Function assumes 2nd column of pd.fn is the celltype phenotype column and renames it "CellType"
# Added a return for probeList to get the top 50 hyper and hypo methylated per cell type (100 in total)
# Returns list: coefEsts = coefEsts, compTable = compTable, sampleMeans = pMeans, probeList = probeList
# Note: this is equivalent to running "both" for probeSelect in original estimeCellCounts minfi implementation

pickCompProbes_kc <- function(p, pd.fn, numProbes = 50) {
  # Custom implementation of pickCompProbes function run with local reference matrices
  # Default for both (50 hyper- and hypo-methylated)

  splitit <- function(x) {
    split(seq_along(x), x)
  }
  
  # Change cell type column label
  colnames(pd.fn)[2] <- "CellType"
  
  ffComp <- rowFtests(p, pd.fn$CellType)
  prof <- vapply(
    X = splitit(pd.fn$CellType),
    FUN = function(j)
      rowMeans2(p, cols = j),
    FUN.VALUE = numeric(nrow(p))
  )
  r <- rowRanges(p)
  compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
  names(compTable)[1] <- "Fstat"
  names(compTable)[c(-2,-1, 0) + ncol(compTable)] <-
    c("low", "high", "range")
  tIndexes <- splitit(pd.fn$CellType)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0, ncol(p))
    x[i] <- 1
    return(rowttests(p, factor(x)))
  })
  
  # Else, in line with non cord blood reference recommendation
  probeList <- lapply(tstatList, function(x) {
    y <- x[x[, "p.value"] < 1e-8,]
    yUp <- y[order(y[, "dm"], decreasing = TRUE),]
    yDown <- y[order(y[, "dm"], decreasing = FALSE),]
    c(rownames(yUp)[seq_len(numProbes)],
      rownames(yDown)[seq_len(numProbes)])
  })
  
  trainingProbes <- unique(unlist(probeList))
  p <- p[trainingProbes, ]
  
  pMeans <- colMeans2(p)
  names(pMeans) <- pd.fn$CellType
  
  form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd.fn$CellType), collapse = "+")))
  phenoDF <- as.data.frame(model.matrix( ~ pd.fn$CellType - 1))
  colnames(phenoDF) <- sub("^pd.fn\\$CellType", "", colnames(phenoDF))
  if (ncol(phenoDF) == 2) {
    # Two group solution
    X <- as.matrix(phenoDF)
    coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
  } else {
    # > 2 groups solution
    tmp <- validationCellType(Y = p,
                              pheno = phenoDF,
                              modelFix = form)
    coefEsts <- tmp$coefEsts
  }
  
  return(list(
    coefEsts = coefEsts,
    compTable = compTable,
    sampleMeans = pMeans,
    probeList = probeList
  ))
}

pickCompProbes_kc_all <- function(p, pd.fn, numProbes = 50) {
  # Custom implementation of pickCompProbes function run with local reference matrices
  # Default for both (50 hyper- and hypo-methylated)
  #numProbes = 50
  
  # Debug
  #p = beta.celltypes
  #pd.fn = pd.celltypes
  #numProbes = dim(beta.celltypes)[1]
  #numProbes = 50
  
  splitit <- function(x) {
    split(seq_along(x), x)
  }
  
  # Change cell type column label
  colnames(pd.fn)[2] <- "CellType"
  
  ffComp <- rowFtests(p, pd.fn$CellType)
  prof <- vapply(
    X = splitit(pd.fn$CellType),
    FUN = function(j)
      rowMeans2(p, cols = j),
    FUN.VALUE = numeric(nrow(p))
  )
  r <- rowRanges(p)
  compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
  names(compTable)[1] <- "Fstat"
  names(compTable)[c(-2,-1, 0) + ncol(compTable)] <-
    c("low", "high", "range")
  tIndexes <- splitit(pd.fn$CellType)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0, ncol(p))
    x[i] <- 1
    return(rowttests(p, factor(x)))
  })
  
  # Else, in line with non cord blood reference recommendation
  probeList <- lapply(tstatList, function(x) {
    y <- x
    yUp <- y[order(y[, "dm"], decreasing = TRUE),]
    yDown <- y[order(y[, "dm"], decreasing = FALSE),]
    c(rownames(yUp)[seq_len(numProbes)],
      rownames(yDown)[seq_len(numProbes)])
  })
  
  trainingProbes <- unique(unlist(probeList))
  p <- p[trainingProbes, ]
  
  pMeans <- colMeans2(p)
  names(pMeans) <- pd.fn$CellType
  
  form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd.fn$CellType), collapse = "+")))
  phenoDF <- as.data.frame(model.matrix( ~ pd.fn$CellType - 1))
  colnames(phenoDF) <- sub("^pd.fn\\$CellType", "", colnames(phenoDF))
  if (ncol(phenoDF) == 2) {
    # Two group solution
    X <- as.matrix(phenoDF)
    coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
  } else {
    # > 2 groups solution
    tmp <- validationCellType(Y = p,
                              pheno = phenoDF,
                              modelFix = form)
    coefEsts <- tmp$coefEsts
  }
  
  return(list(
    coefEsts = coefEsts,
    compTable = compTable,
    sampleMeans = pMeans,
    probeList = probeList
  ))
}

# Custom implementation of estimateLC function run ewastools deconvolution with local reference matrices
estimateLC_kc = function(meth, ref, constrained=TRUE){
  
  J = ncol(meth)
  
  coefs = read.table(ref)
  coefs = as.matrix(coefs)
  n_celltypes = ncol(coefs)
  
  markers = match(rownames(coefs),rownames(meth))
  EST = sapply(1:J,function(j){
    tmp = meth[markers,j]
    i = !is.na(tmp)
    
    if(constrained == FALSE){
      return(
        quadprog::solve.QP(
          t(coefs[i,]) %*% coefs[i,]
          ,t(coefs[i,]) %*% tmp[i]
          ,diag(n_celltypes)
          ,rep(0,n_celltypes)
        )$sol)
    }else{
      return(
        quadprog::solve.QP(
          t(coefs[i,]) %*% coefs[i,]
          ,t(coefs[i,]) %*% tmp[i]
          ,cbind(1,diag(n_celltypes))
          ,c(1,rep(0,n_celltypes))
          ,meq=1
        )$sol)
    }
  })
  EST = t(EST)
  colnames(EST) = colnames(coefs)
  EST = data.table(EST)
  
  return(EST)
}

# Purely for debuggin purposes
pickCompProbes_kc_debug <- function(p, pd.fn, numProbes = 50) {
  # Custom implementation of pickCompProbes function run with local reference matrices
  # Default for both (50 hyper- and hypo-methylated)
  #numProbes = 50
  
  p = beta.celltypes
  pd.fn = pd.celltypes
  numProbes = 50
  
  splitit <- function(x) {
    split(seq_along(x), x)
  }
  
  # Change cell type column label
  colnames(pd.fn)[2] <- "CellType"
  
  ffComp <- rowFtests(p, pd.fn$CellType)
  prof <- vapply(
    X = splitit(pd.fn$CellType),
    FUN = function(j)
      rowMeans2(p, cols = j),
    FUN.VALUE = numeric(nrow(p))
  )
  r <- rowRanges(p)
  compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
  names(compTable)[1] <- "Fstat"
  names(compTable)[c(-2,-1, 0) + ncol(compTable)] <-
    c("low", "high", "range")
  tIndexes <- splitit(pd.fn$CellType)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0, ncol(p))
    x[i] <- 1
    return(rowttests(p, factor(x)))
  })
  
  # Else, in line with non cord blood reference recommendation
  probeList <- lapply(tstatList, function(x) {
    y <- x[x[, "p.value"] < 1e-8,]
    yUp <- y[order(y[, "dm"], decreasing = TRUE),]
    yDown <- y[order(y[, "dm"], decreasing = FALSE),]
    c(rownames(yUp)[seq_len(numProbes)],
      rownames(yDown)[seq_len(numProbes)])
  })
  
  trainingProbes <- unique(unlist(probeList))
  p <- p[trainingProbes, ]
  
  pMeans <- colMeans2(p)
  names(pMeans) <- pd.fn$CellType
  
  form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd.fn$CellType), collapse = "+")))
  phenoDF <- as.data.frame(model.matrix( ~ pd.fn$CellType - 1))
  colnames(phenoDF) <- sub("^pd.fn\\$CellType", "", colnames(phenoDF))
  if (ncol(phenoDF) == 2) {
    # Two group solution
    X <- as.matrix(phenoDF)
    coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
  } else {
    # > 2 groups solution
    tmp <- validationCellType(Y = p,
                              pheno = phenoDF,
                              modelFix = form)
    coefEsts <- tmp$coefEsts
  }
  
  return(list(
    coefEsts = coefEsts,
    compTable = compTable,
    sampleMeans = pMeans
  ))
}

# Custom implementation of estimateLC function run ewastools deconvolution with local reference matrices
estimateLC_kc = function(meth, ref, constrained=TRUE){
  
  J = ncol(meth)
  
  coefs = read.table(ref)
  coefs = as.matrix(coefs)
  n_celltypes = ncol(coefs)
  
  markers = match(rownames(coefs),rownames(meth))
  EST = sapply(1:J,function(j){
    tmp = meth[markers,j]
    i = !is.na(tmp)
    
    if(constrained == FALSE){
      return(
        quadprog::solve.QP(
          t(coefs[i,]) %*% coefs[i,]
          ,t(coefs[i,]) %*% tmp[i]
          ,diag(n_celltypes)
          ,rep(0,n_celltypes)
        )$sol)
    }else{
      return(
        quadprog::solve.QP(
          t(coefs[i,]) %*% coefs[i,]
          ,t(coefs[i,]) %*% tmp[i]
          ,cbind(1,diag(n_celltypes))
          ,c(1,rep(0,n_celltypes))
          ,meq=1
        )$sol)
    }
  })
  EST = t(EST)
  colnames(EST) = colnames(coefs)
  EST = data.table(EST)
  
  return(EST)
}


save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

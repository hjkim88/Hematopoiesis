###
#   File name : Additional_Analysis_Nov_2020.R
#   Author    : Hyunjin Kim
#   Date      : Nov 18, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Perform additional analyses for Trent's project based on the
#               "20201119 RNA-Magnet notes for Hyunjin.txt".
#   
#   Instruction
#               1. Source("Additional_Analysis_Oct_2020.R")
#               2. Run the function "additional_analysis" - specify the input paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Additional_Analysis_Oct_2020.R/Additional_Analysis_Oct_2020.R")
#               > additional_analysis(Robj_path="./data/Combined_Seurat_Obj.RDS",
#                                     outputDir="./results/Additional_Nov2020/")
###

additional_analysis <- function(Robj_path="./data/Combined_Seurat_Obj.RDS",
                                outputDir="./results/Additional_Nov2020/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(ggrepel, quietly = TRUE)) {
    install.packages("ggrepel")
    require(ggrepel, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(RNAMagnet, quietly = TRUE)) {
    remotes::install_github("veltenlab/rnamagnet")
    require(RNAMagnet, quietly = TRUE)
  }
  # Sys.setenv(RETICULATE_MINICONDA_PATH = "D:/Miniconda3/")
  # Sys.setenv(RETICULATE_PYTHON = "D:/Miniconda3/envs/r-reticulate/")
  if(!require(reticulate, quietly = TRUE)) {
    install.packages("reticulate")
    require(reticulate, quietly = TRUE)
  }
  if(!require(Rmagic, quietly = TRUE)) {
    install.packages("Rmagic")
    require(Rmagic, quietly = TRUE)
  }
  if(!require(scales, quietly = TRUE)) {
    install.packages("scales")
    library(scales, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  if(!require(ggbeeswarm, quietly = TRUE)) {
    install.packages("ggbeeswarm")
    require(ggbeeswarm, quietly = TRUE)
  }
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
  }
  if(!require(RColorBrewer, quietly = TRUE)) {
    install.packages("RColorBrewer")
    require(RColorBrewer, quietly = TRUE)
  }
  if(!require(gplots, quietly = TRUE)) {
    install.packages("gplots")
    library(gplots, quietly = TRUE)
  }
  
  ### result directory
  outputDir2 <- paste0(outputDir, "RNA_Magnet/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### In terminal, if 'bash: conda: command not found',
  ### . /d/Miniconda3/etc/profile.d/conda.sh
  ###
  ### python --version
  ### python-magic needs python version 2.7 not 3.*
  ### conda search python
  ### conda install python=2.7
  ###
  ### To remove an environment
  ### conda deactivate
  ### conda remove --name "r-reticulate" --all
  ###
  ### conda create -n "r-reticulate"
  ### conda install -n "r-reticulate" python=2.7
  ### conda install -n "r-reticulate" numpy=1.15.0
  ###
  ### conda activate "r-reticulate"
  ### pip install python-magic=0.4.18
  ### conda install -c conda-forge python-magic
  ###
  ### pip list
  ###
  ### https://github.com/pidydx/libmagicwin64
  ### Drop the dlls (libgnurx-0.dll, magic.mgc, & magic1.dll) in C:\Windows\System32 and python magic will import correctly.
  
  ### see python environment since RNAMagnet uses the python module 'magic'
  # conda search python-magic --channel conda-forge
  # conda_create(packages = c("python=3.8.5", "numpy=1.19.1"))
  # conda_install("r-reticulate", "magic-impute")
  # pip install magic-impute
  # there are two magics: magic-impute & python-magic.
  # python-magic is not what we want, so don't install it - magic-impute would work fine.
  use_condaenv("r-reticulate")
  py_config()
  py_module_available("magic")
  
  ### updated stroma RDS file
  Updated_Seurat_Obj <- readRDS(file = Robj_path)
  
  ### create new column for the analysis
  Updated_Seurat_Obj@meta.data$Dev_Anno <- paste0(Updated_Seurat_Obj@meta.data$Development, "_",
                                                  Updated_Seurat_Obj@meta.data$Annotation)
  
  ### determine necessary variables
  time_points <- c("E16", "E18", "P0", "ADULT")
  HSPC_populations <- c("LTHSC", "STHSC", "MPP2", "MPP3", "MPP4")
  
  ### set the ident of the object with Tissue type
  Updated_Seurat_Obj <- SetIdent(object = Updated_Seurat_Obj,
                                 cells = rownames(Updated_Seurat_Obj@meta.data),
                                 value = Updated_Seurat_Obj@meta.data$Tissue)
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Updated_Seurat_Obj)), rownames(Updated_Seurat_Obj@meta.data)))
  
  ### run UMAP
  Updated_Seurat_Obj <- FindVariableFeatures(Updated_Seurat_Obj)
  Updated_Seurat_Obj <- ScaleData(Updated_Seurat_Obj)
  Updated_Seurat_Obj <- RunUMAP(Updated_Seurat_Obj, dims = 1:5)
  
  ### a function for color brewer
  cell_pal <- function(cell_vars, pal_fun) {
    if (is.numeric(cell_vars)) {
      pal <- pal_fun(100)
      return(pal[cut(cell_vars, breaks = 100)])
    } else {
      categories <- sort(unique(cell_vars))
      pal <- setNames(pal_fun(length(categories)), categories)
      return(pal[cell_vars])
    }
  }
  
  ### check the UMAP for all the cells
  umap_plot <- DimPlot(Updated_Seurat_Obj, reduction = "umap", group.by = "HSPC", pt.size = 1.5) +
    labs(title = paste0("UMAP_ALL_CELL_TYPES"))
  umap_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  plot(umap_plot)
  
  umap_plot <- DimPlot(Updated_Seurat_Obj, reduction = "umap", group.by = "Cell_Type", pt.size = 1.5) +
    labs(title = paste0("UMAP_ALL_CELL_TYPES"))
  umap_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  plot(umap_plot)
  
  ### set the ident of the object with the group info
  Updated_Seurat_Obj <- SetIdent(object = Updated_Seurat_Obj,
                                 cells = rownames(Updated_Seurat_Obj@meta.data),
                                 value = Updated_Seurat_Obj@meta.data$Development)
  
  ### split the Seurat obj based on the given info
  Combined_Adult_Seurat_Obj <- subset(Updated_Seurat_Obj, idents="ADULT")
  
  ### draw a UMAP with cell type
  umap_plot <- DimPlot(Combined_Adult_Seurat_Obj, reduction = "umap", group.by = "Cell_Type", pt.size = 1.5) +
    labs(title = paste0("UMAP_Combined_Adults_Only"))
  umap_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  plot(umap_plot)
  ggsave(file = paste0(outputDir2, "UMAP_Combined_Adult_Cells.png"), width = 15, height = 10, dpi = 300)
  
  ### draw a UMAP with HSPC info
  umap_plot <- DimPlot(Combined_Adult_Seurat_Obj, reduction = "umap", group.by = "HSPC", pt.size = 1.5) +
    labs(title = paste0("UMAP_Combined_Adults_Only"))
  umap_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  plot(umap_plot)
  ggsave(file = paste0(outputDir2, "UMAP_Combined_Adult_HSPC.png"), width = 15, height = 10, dpi = 300)
  
  ### draw a UMAP with Trent's annotation
  Combined_Adult_Seurat_Obj@meta.data$Annotation <- factor(Combined_Adult_Seurat_Obj@meta.data$Annotation)
  umap_plot <- DimPlot(Combined_Adult_Seurat_Obj, reduction = "umap", group.by = "Annotation", pt.size = 1.5) +
    labs(title = paste0("UMAP_Combined_Adult_With_the_Annotation"))
  umap_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  umap_plot <- LabelClusters(plot = umap_plot, id = "Annotation", col = "black")
  plot(umap_plot)
  ggsave(file = paste0(outputDir2, "UMAP_Combined_Adult_With_the_Annotation.png"), width = 15, height = 10, dpi = 300)
  
  
  #'Low-level function to run core RNA magnet steps
  #'
  #' Users are advised to use the top-level functions \code{\link{RNAMagnetAnchors}} and \code{\link{RNAMagnetSignaling}} which appropriately set default parameters and return user-friendly return values.
  #' This is a low level function for development purposes.
  #'@param seurat An object of class \code{\link[Seurat]{seurat}} containing a valid clustering and t-SNE information
  #'
  #'. For information on how to create such an object, see https://satijalab.org/seurat/get_started.html
  #'@param anchors A character vector of anchor populations. Entries must be levels of the seurat identity vector. If \code{NULL}: All entries of the seurat identity vector are used as anchors.
  #'@param neighborhood.distance See \code{\link{RNAMagnetAnchors}}
  #'@param neighborhood.gradient See \code{\link{RNAMagnetAnchors}}
  #'@param .k Fuzzification parameter, see detail. Recommended to leave at the default value.
  #'@param .x0 Fuzzification parameter, see detail. Recommended to leave at the default value.
  #'@param .minExpression Minimal expression level of genes to be included, specified as the number of cells in the dataset that express the gene.
  #'@param .minMolecules Number of molecules per cell required to count the cells as positive.
  #'@param .version The version of the underlying ligand-receptor database. See \code{\link{getLigandsReceptors}}.
  #'@param .cellularCompartment Types of ligands to be included. For physical interactions, defaults to \code{ c("Membrane","ECM","Both")}.  See \code{\link{getLigandsReceptors}}.
  #'@param .manualAnnotation Annotation status of ligands to be included. Default to \code{"Correct"}. See \code{\link{getLigandsReceptors}}.
  #'@param .symmetric Assume that if A is a receptor for B, B is also a receptor for A
  #'@details The algorithm takes the following steps: \enumerate{
  #'\item Ligand-receptor pairs are selected based on the parameters \code{.version}, \code{.cellularCompartment} and \code{.manualAnnotation}. Choice of \code{.cellularCompartment} is crucial for determining the algorithm's behavior, e.g. if set to \code{c("Secreted","Both")}, paracrine signaling interactions involving soluble ligands are investigated.
  #'\item Dropout values in the expression levels of ligands and receptors are imputed using \code{\link[Rmagic]{magic}}
  #'\item Mean expression level of ligands and receptors is computed for all anchor populations
  #'\item For each cell or anchor population, the expression of each ligand and receptor is encoded as a fuzzy logic variable
  #'\item Fuzzy logic AND is used to compute probabilities for a given interaction to be active between a single cell and an anchor population
  #'\item An interaction score is computed as the sum of interaction probabilities across all possible ligand-receptor pairs
  #'\item Specificty scores are computed by comparing interaction scores to average interaction scores in a local neighborhood.
  #'}
  #'@details Add the methods section of the paper here!
  #'@return Returns an object of class \code{\link{rnamagnet}}
  #'@export
  RNAMagnetBase <- function(seurat, anchors=NULL,neighborhood.distance=NULL, neighborhood.gradient =NULL, .k = 10, .x0 = 0.5, .minExpression,.minMolecules=1, .version = "1.0.0", .cellularCompartment, .manualAnnotation = "Correct", .symmetric = F) {
    cat("Setting everything up...\n")
    
    if (grepl("^3", Biobase::package.version("Seurat"))) {
      seurat.ident <- Idents(seurat)
      seurat.cell.names <- colnames(seurat)
      seurat.pca <- Embeddings(seurat, reduction = "pca")
      seurat.raw.data <- GetAssayData(seurat, slot = "counts")
      
    } else {
      seurat.ident <- seurat@ident
      seurat.cell.names <- seurat@cell.names
      seurat.pca <- seurat@dr$pca@cell.embeddings
      seurat.raw.data <- seurat@raw.data
    }
    
    if (is.null(anchors)) anchors <- as.character(unique(seurat.ident))
    
    out <- new("rnamagnet", celltype = seurat.ident, params = list("neighborhood.distance"=neighborhood.distance, "neighborhood.gradient" =neighborhood.gradient, ".k" = .k, ".x0" = .x0, ".minExpression" = .minExpression, ".minMolecules" = .minMolecules, ".cellularCompartment" = .cellularCompartment, ".manualAnnotation" = .manualAnnotation, ".symmetric" = .symmetric))
    
    #compute cell-cell similarity
    similarity <- as.matrix(1-cor(t(seurat.pca[,1:15])))
    
    #prepare database
    ligrec <- getLigandsReceptors(.version, .cellularCompartment, .manualAnnotation)
    if (.symmetric) ligrec <- makeSymmetric(ligrec) #for physical interactions: i a binds b, b binds a.
    
    
    #prepare genes included into MAGIC
    filteredGenes <- rownames(seurat.raw.data)[apply(seurat.raw.data[,seurat.cell.names] >=.minMolecules ,1,sum) > .minExpression]
    
    genes <- unique(c(ligrec$Receptor.Mouse,
                      ligrec$Ligand.Mouse))
    
    genes <- genes[sapply(genes, function(x) {
      entries <- strsplit(x, "[&|]")
      if (grepl("&", x))
        all(entries[[1]] %in% filteredGenes)
      else any(entries[[1]] %in% filteredGenes)
    })]
    
    genes_formagic <- unlist(strsplit( genes,"[&|]"))
    
    
    formagic <- Matrix::t(seurat.raw.data[,seurat.cell.names])
    formagic <- Rmagic::library.size.normalize(formagic)
    formagic <- sqrt(formagic)
    
    #Run MAGIC
    cat("Now running MAGIC to impute dropout values...\n")
    mymagic <- Rmagic::magic(formagic, genes = genes_formagic, seed = 0xbeef)
    mymagic <- as.data.frame(mymagic)
    
    #handle expression levels of heterodimers
    resolvedRawData <- resolveData(t(mymagic), genes)
    resolvedRawData <- t(apply(resolvedRawData, 1, function(x) (x - min(x )) / (max(x) - min(x)) ))
    
    
    annotated_genes <- rownames(resolvedRawData)
    out@mylr <- subset(ligrec, Receptor.Mouse %in% annotated_genes &
                         Ligand.Mouse %in% annotated_genes)
    
    stepf<- Vectorize(function(x) if (x<0) 0 else x)
    
    cat("Now running RNAMagnet...\n")
    
    #compute mean gene expression level per population
    out@anchors <- do.call(cbind, lapply(anchors, function(id) {
      apply(resolvedRawData[,seurat.ident == id,drop=FALSE],1,mean)
    }));
    colnames(out@anchors) <- anchors
    
    #use fuzzy logic and operations to compute interaction score
    out@interaction <- sapply(anchors, function(pop_l) {
      out@mylr$expression_ligand <- out@anchors[out@mylr$Ligand.Mouse, pop_l]
      sapply(seurat.cell.names, function(cell_r) {
        out@mylr$expression_receptor <-resolvedRawData[out@mylr$Receptor.Mouse, cell_r]
        sum(kernel(out@mylr$expression_ligand, k =.k, x0 = .x0) * kernel(out@mylr$expression_receptor, k = .k, x0=.x0 )) #kernel performs fuzzification
      })
    })
    ### sometimes this is missing
    colnames(out@interaction) <- anchors
    
    out@specificity <- t(sapply(rownames(out@interaction), function(cell) {
      x <- out@interaction[cell,]
      beta <- x/sum(x)
      if (!is.null(neighborhood.distance)) alpha <- apply(out@interaction * (1-kernel(similarity[cell,],neighborhood.gradient,x0=neighborhood.distance )),2,sum) else alpha <- apply(out@interaction,2,mean)
      alpha <- alpha / sum(alpha)
      beta- alpha
    }))
    rownames(out@specificity) <- rownames(out@interaction)
    
    out@adhesiveness <- apply(mymagic, 1,function(x) sum(kernel(x,x0 = .x0, k = .k)))
    return(out)
    
  }
  
  ### RNAMagnetAnchors
  RNAMagnetAnchors <- function(seurat, anchors, return = "summary", neighborhood.distance = 0.7, neighborhood.gradient = 3, .k = 10, .x0 = 0.5, .minExpression = 0, .minMolecules = 1, .version = "1.0.0", .cellularCompartment = c("Membrane","ECM","Both"), .manualAnnotation = "Correct" ) {
    
    myMagnet <- RNAMagnetBase(seurat, anchors, neighborhood.distance,neighborhood.gradient, .k, .x0, .minExpression, .minMolecules, .version, .cellularCompartment, .manualAnnotation,TRUE)
    if (return=="rnamagnet-class") myMagnet else data.frame(direction = as.factor(colnames(myMagnet@specificity)[apply(myMagnet@specificity,1,which.max)]), adhesiveness = myMagnet@adhesiveness, myMagnet@specificity[,anchors])
    
  }
  
  ### hierarchical clustering functions
  dist.spear <- function(x) as.dist(1-cor(t(x), method = "spearman"))
  hclust.ave <- function(x) hclust(x, method="average")
  
  ### A function that creates mupliple plots with highliting
  ### one anchor at a time - thus result in many plots
  ### only one anchor is colored and others will be grayed
  ### + heatmap that Jeremy suggested
  ### plot_df example below:
  # plot_df <- data.frame(X=dim_map[rownames(SO@meta.data),1],
  #                       Y=dim_map[rownames(SO@meta.data),2],
  #                       cluster_color = SO@meta.data[,target_col],
  #                       direction = SO@meta.data$direction,
  #                       direction2 = SO@meta.data$direction2,
  #                       adhesiveness = SO@meta.data$adhesiveness,
  #                       adhesiveness2 = SO@meta.data$adhesiveness2,
  #                       specificity = SO@meta.data$specificity,
  #                       specificity2 = SO@meta.data$specificity2,
  #                       stringsAsFactors = FALSE, check.names = FALSE)
  ### original_col: the column name of plot_df that shows the original cell type
  ### anchor_col: the column name of plot_df that will be colored in the plots (Direction)
  ### transp_col: adhesiveness/specificity column of the plot_df that will be represented with transparency
  ### dim_method: Is the RNAMagent based on UMAP or PCA? Should be one of them
  ### title: title of the plots
  ### output_directory: the directory that the results of the function will be saved
  ### width, height, dpi: plot width, height, and resolution values
  ###
  ### Example
  ### original_col = "cluster_color2"
  ### anchor_col = "direction2"
  ### transp_col = "adhesiveness2"
  RNAMagnet_multiple_anchor_plot <- function(plot_df,
                                             original_col,
                                             anchor_col,
                                             transp_col,
                                             dim_method = c("UMAP", "PCA"),
                                             title,
                                             output_directory = "./") {
    
    ### get colors for the clustering result
    original_colors_clust <- cell_pal(unique(as.character(plot_df[,original_col])), hue_pal())
    anchor_colors_clust <- original_colors_clust[unique(as.character(plot_df[,anchor_col]))]
    
    ### PCA/UMAP label
    if(dim_method[1] == "PCA") {
      x_lab <- "PC1"
      y_lab <- "PC2"
    } else if(dim_method[1] == "UMAP") {
      x_lab <- "UMAP1"
      y_lab <- "UMAP2"
    } else {
      stop("ERROR: dim_method not PCA nor UMAP.")
    }
    
    ### Adhresiveness/Specificity label
    if(grepl("adhesiveness", transp_col)) {
      transp_label <- "Adhesiveness"
    } else if(grepl("specificity", transp_col)) {
      transp_label <- "Specificity"
    } else {
      stop("ERROR: transp_col not adhesiveness nor specificity")
    }
    
    ### make the columns as factors
    plot_df[,original_col] <- factor(plot_df[,original_col], levels = unique(plot_df[,original_col]))
    plot_df[,anchor_col] <- factor(plot_df[,anchor_col], levels = unique(plot_df[,anchor_col]))
    
    ### set this globally for LabelClusters
    options(ggrepel.max.overlaps = Inf)
    
    ### create an empty list
    p <- vector("list", length = length(unique(plot_df[,anchor_col]))+2)
    
    #
    ### draw RNAMagnet plots coloring differently for each anchor
    #
    ### original
    p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col=original_col), size=2) +
      ggtitle(paste0("Original ", dim_method[1])) +
      xlab(x_lab) + ylab(x_lab) +
      labs(col="Cell Type") +
      theme_classic(base_size = 16) +
      scale_color_manual(values = original_colors_clust,
                         labels = names(original_colors_clust)) +
      theme(legend.title = element_text(size = 10),
            legend.text = element_text(size = 8))
    ### the id column in the plot_df should be a factor
    p[[1]] <- LabelClusters(plot = p[[1]], id = original_col, col = "black")
    
    ### with all the directions
    p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col=anchor_col, alpha=transp_col), size=2) +
      ggtitle(paste0("RNAMagnet with all the directions")) +
      xlab(x_lab) + ylab(x_lab) +
      labs(col="Direction", alpha=transp_label) +
      theme_classic(base_size = 16) +
      scale_color_manual(values = anchor_colors_clust,
                         labels = names(anchor_colors_clust)) +
      theme(legend.title = element_text(size = 10),
            legend.text = element_text(size = 8))
    ### the id column in the plot_df should be a factor
    p[[2]] <- LabelClusters(plot = p[[2]], id = original_col, col = "black")
    for(i in 1:length(unique(plot_df[,anchor_col]))) {
      specific_colorset <- anchor_colors_clust[as.character(plot_df[,anchor_col])]
      specific_colorset[which(specific_colorset != anchor_colors_clust[i])] <- "gray"
      plot_df$specific_colorset <- specific_colorset
      p[[i+2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
        geom_point(aes_string(col="specific_colorset", alpha=transp_col), size=2) +
        ggtitle(paste0("RNAMagnet (Coloring ", names(anchor_colors_clust)[i], " only)")) +
        xlab(x_lab) + ylab(x_lab) +
        labs(col="Direction", alpha=transp_label) +
        theme_classic(base_size = 16) +
        scale_color_manual(values = c(anchor_colors_clust[i], "gray"),
                           labels = c(names(anchor_colors_clust)[i], "Others")) +
        theme(legend.title = element_text(size = 10),
              legend.text = element_text(size = 8))
      ### the id column in the plot_df should be a factor
      p[[i+2]] <- LabelClusters(plot = p[[i+2]], id = original_col, col = "black")
    }
    
    ### save the plots
    if(length(p) > 20 && length(p) <= 25) {
      fig_nrow <- 5
      fig_ncol <- 5
    } else if(length(p) > 16 && length(p) <= 20) {
      fig_nrow <- 4
      fig_ncol <- 5
    } else if(length(p) > 12 && length(p) <= 16) {
      fig_nrow <- 4
      fig_ncol <- 4
    } else if(length(p) > 9 && length(p) <= 12) {
      fig_nrow <- 3
      fig_ncol <- 4
    } else if(length(p) > 6 && length(p) <= 9) {
      fig_nrow <- 3
      fig_ncol <- 3
    } else if(length(p) > 4 && length(p) <= 6) {
      fig_nrow <- 2
      fig_ncol <- 3
    } else if(length(p) > 2 && length(p) <= 4) {
      fig_nrow <- 2
      fig_ncol <- 2
    } else if(length(p) > 1 && length(p) <= 2) {
      fig_nrow <- 1
      fig_ncol <- 2
    } else if(length(p) == 1) {
      fig_nrow <- 1
      fig_ncol <- 1
    } else {
      fig_nrow <- 6
      fig_ncol <- 6
    }
    g <- arrangeGrob(grobs = p,
                     nrow = fig_nrow,
                     ncol = fig_ncol,
                     top = paste0(dim_method[1], " with Direction & ", transp_label))
    ggsave(file = paste0(output_directory, title, "_", dim_method[1], "_", anchor_col, "_and_", transp_col, "_Scatter.png"),
           g, width = 20, height = 12, dpi = 300)
    
    
    ### heatmap that Jeremy suggested
    ### rows: original cell type
    ### columns: direction
    ### color: direction
    ### transparency: adhesiveness/specificity
    #
    ### make a matrix for the heatmap
    heatmap_mat <- matrix(0, nrow(plot_df), length(unique(as.character(plot_df[,original_col]))))
    rownames(heatmap_mat) <- rownames(plot_df)
    colnames(heatmap_mat) <- unique(as.character(plot_df[,original_col]))
    
    ### fill the matrix
    for(i in 1:nrow(heatmap_mat)) {
      heatmap_mat[i,as.character(plot_df[i,anchor_col])] <- as.numeric(plot_df[i,transp_col])
    }
    
    ### set rowside colors
    uniqueV <- unique(as.character(plot_df[,original_col]))
    colors <- colorRampPalette(brewer.pal(length(uniqueV), "Spectral"))(length(uniqueV))
    names(colors) <- uniqueV
    
    ### heatmap
    png(paste0(output_directory, title, "_", dim_method[1], "_", anchor_col, "_and_", transp_col, "_Heatmap.png"),
        width = 2000, height = 1200, res = 200)
    par(oma=c(3,0,3,0))
    heatmap.2(heatmap_mat,
              xlab = "", ylab = "", col=colorpanel(10, low = "white", high = "black"),
              scale="none", key=TRUE, key.title = transp_label, keysize=1.5, density.info="density",
              dendrogram = "none", trace = "none",
              labRow = FALSE, labCol = colnames(heatmap_mat),
              Rowv = TRUE, Colv = FALSE,
              distfun=dist.spear, hclustfun=hclust.ave,
              RowSideColors = rbind(colors[as.character(plot_df[,original_col])]),
              ColSideColors = cbind(colors[as.character(colnames(heatmap_mat))]),
              cexRow = 1.5, cexCol = 1.5, na.rm = TRUE,
              main = paste0(title, "_(",
                            nrow(heatmap_mat), " Cells x ",
                            ncol(heatmap_mat), " Directions)", "\n",
                            "[Rows: Original Cell Type, Columns: Direction]"))
    legend("left", inset = 0, xpd = TRUE, title = "Cell Type", legend = names(colors), fill = colors, cex = 1, box.lty = 0)
    dev.off()
    
  }
  
  ### RNAMagnet for:
  ### Adult LT-HSCs vs Adult CARs
  ### Adult LT-HSCs vs Adult SECs
  ### Adult LT-HSCs vs Adult AECs
  ### Adult LT-HSCs vs Adult F1
  ### Adult LT-HSCs vs Adult F4
  ### Adult LT-HSCs vs Adult F5
  
  ### put two annotation types into one
  Combined_Adult_Seurat_Obj@meta.data$Annotation2 <- Combined_Adult_Seurat_Obj@meta.data$Cell_Type
  Combined_Adult_Seurat_Obj@meta.data$Annotation2[which(Combined_Adult_Seurat_Obj@meta.data$HSPC == "LTHSC")] <- "LTHSC"
  Combined_Adult_Seurat_Obj@meta.data$Annotation2[which(Combined_Adult_Seurat_Obj@meta.data$Annotation2 == "Stroma")] <- as.character(Combined_Adult_Seurat_Obj@meta.data$Annotation[which(Combined_Adult_Seurat_Obj@meta.data$Annotation2 == "Stroma")])
  
  ### simply run RNAMagnet only to get signaling interactions
  ### SO: Seurat object
  ### target_col: Target column name of the given seurat object's meta.data
  ### conds: a character vector of all the conditions in the 'target_col' that will be used
  ### result_dir: output directory
  simple_RNAMagnet_interaction <- function(SO, target_col, conds, dim_method, result_dir) {
    
    result_dir <- paste0(result_dir, "/", paste(conds, collapse = "_vs_"), "/")
    dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
    
    ### set the ident of the object with the HSPC type
    SO <- SetIdent(object = SO,
                   cells = rownames(SO@meta.data),
                   value = SO@meta.data[,target_col])
    
    ### split the Seurat obj based on the given info
    SO <- subset(SO, idents=conds)
    
    ### run PCA/UMAP
    if(dim_method == "PCA") {
      SO <- RunPCA(SO, npcs = 15)
      dim_map <- Embeddings(SO, reduction = "pca")[rownames(SO@meta.data),]
      x_lab <- "PC1"
      y_lab <- "PC2"
    } else if(dim_method == "UMAP") {
      SO <- RunUMAP(SO, dims = 1:15)
      dim_map <- Embeddings(SO, reduction = "umap")[rownames(SO@meta.data),]
      x_lab <- "UMAP1"
      y_lab <- "UMAP2"
    } else {
      stop("ERROR: dim_method not PCA nor UMAP.")
    }
    
    ### run RNAMagnet with anchors
    result <- RNAMagnetAnchors(SO,
                               anchors = unique(SO@meta.data[,target_col]))
    
    ### run RNAMagnet with different anchors (new annotation from Trent)
    SO@meta.data$Annotation <- factor(SO@meta.data$Annotation, levels = unique(SO@meta.data$Annotation))
    SO <- SetIdent(object = SO,
                   cells = rownames(SO@meta.data),
                   value = SO@meta.data$Annotation)
    result2 <- RNAMagnetAnchors(SO, anchors = unique(SO@meta.data$Annotation))
    
    ### write the result as an Excel file
    write.xlsx2(data.frame(Cell=rownames(result), result,
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(result_dir, paste(conds, collapse = "_vs_"),
                              "_RNAMagnet_Result.xlsx"),
                sheetName = paste0("RNAMagnet_Result"),
                row.names = FALSE)
    write.xlsx2(data.frame(Cell=rownames(result2), result2,
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(result_dir, paste(conds, collapse = "_vs_"),
                              "_RNAMagnet_Result_Using_Annotation.xlsx"),
                sheetName = paste0("RNAMagnet_Result"),
                row.names = FALSE)
    
    ### add RNAMagnet info to the seurat object
    SO@meta.data$direction <- as.character(result[rownames(SO@meta.data),"direction"])
    SO@meta.data$direction2 <- as.character(result2[rownames(SO@meta.data),"direction"])
    SO@meta.data$adhesiveness <- as.numeric(result[rownames(SO@meta.data),"adhesiveness"])
    SO@meta.data$adhesiveness2 <- as.numeric(result2[rownames(SO@meta.data),"adhesiveness"])
    SO@meta.data$specificity <- as.numeric(sapply(rownames(SO@meta.data), function(x) {
      result[x, result[x,"direction"]]
    }))
    SO@meta.data$specificity2 <- as.numeric(sapply(rownames(SO@meta.data), function(x) {
      result2[x, result2[x,"direction"]]
    }))
    
    ### make a data frame for ggplot
    plot_df <- data.frame(X=dim_map[rownames(SO@meta.data),1],
                          Y=dim_map[rownames(SO@meta.data),2],
                          cluster_color = SO@meta.data[,target_col],
                          cluster_color2 = SO@meta.data$Annotation,
                          direction = SO@meta.data$direction,
                          direction2 = SO@meta.data$direction2,
                          adhesiveness = SO@meta.data$adhesiveness,
                          adhesiveness2 = SO@meta.data$adhesiveness2,
                          specificity = SO@meta.data$specificity,
                          specificity2 = SO@meta.data$specificity2,
                          stringsAsFactors = FALSE, check.names = FALSE)
    
    ### make a scatter plot and a heatmap with the 'plot_df'
    #
    ### Adhesiveness with the original aanotation
    RNAMagnet_multiple_anchor_plot(plot_df = plot_df,
                                   original_col = "cluster_color",
                                   anchor_col = "direction",
                                   transp_col = "adhesiveness",
                                   dim_method = "UMAP",
                                   title = paste0(paste(conds, collapse = "_vs_"),
                                                  "_RNAMagnet_AD"),
                                   output_directory = result_dir)
    ### Specificity with the original aanotation
    RNAMagnet_multiple_anchor_plot(plot_df = plot_df,
                                   original_col = "cluster_color",
                                   anchor_col = "direction",
                                   transp_col = "specificity",
                                   dim_method = "UMAP",
                                   title = paste0(paste(conds, collapse = "_vs_"),
                                                  "_RNAMagnet_SP"),
                                   output_directory = result_dir)
    ### Adhesiveness with Trent's aanotation
    RNAMagnet_multiple_anchor_plot(plot_df = plot_df,
                                   original_col = "cluster_color2",
                                   anchor_col = "direction2",
                                   transp_col = "adhesiveness2",
                                   dim_method = "UMAP",
                                   title = paste0(paste(conds, collapse = "_vs_"),
                                                  "_RNAMagnet_AD"),
                                   output_directory = result_dir)
    ### Specificity with Trent's aanotation
    RNAMagnet_multiple_anchor_plot(plot_df = plot_df,
                                   original_col = "cluster_color2",
                                   anchor_col = "direction2",
                                   transp_col = "specificity2",
                                   dim_method = "UMAP",
                                   title = paste0(paste(conds, collapse = "_vs_"),
                                                  "_RNAMagnet_SP"),
                                   output_directory = result_dir)
    
    # ### scatter plot
    # p <- list()
    # 
    # ### draw a scatter plot with the adhesiveness info
    # p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    #   geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
    #   xlab(x_lab) + ylab(x_lab) +
    #   labs(col="Cluster") +
    #   ggtitle(paste0("UMAP with Cell Type")) +
    #   scale_color_brewer(palette="Dark2") +
    #   theme_classic(base_size = 16)
    # 
    # ### get colors for the clustering result
    # cell_colors_clust <- cell_pal(unique(SO@meta.data[,target_col]), hue_pal())
    # 
    # p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    #   geom_point(aes_string(col="direction", alpha="adhesiveness"), size=2) +
    #   xlab(x_lab) + ylab(x_lab) +
    #   labs(col="Direction", alpha="Adhesiveness") +
    #   ggtitle(paste("UMAP with Direction & Adhesiveness")) +
    #   theme_classic(base_size = 16) +
    #   scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$direction)[order(unique(plot_df$direction))])],
    #                      labels = names(cell_colors_clust[as.character(unique(plot_df$direction)[order(unique(plot_df$direction))])]))
    # 
    # ### get colors for the clustering result
    # cell_colors_clust <- cell_pal(unique(SO@meta.data$Annotation), hue_pal())
    # 
    # plot_df$direction2 <- factor(plot_df$direction2, levels = unique(plot_df$direction2))
    # p[[3]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    #   geom_point(aes_string(col="direction2", alpha="adhesiveness2"), size=2) +
    #   xlab(x_lab) + ylab(x_lab) +
    #   labs(col="Direction", alpha="Adhesiveness") +
    #   ggtitle("UMAP with Direction & Adhesiveness") +
    #   theme_classic(base_size = 16) +
    #   scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])],
    #                      labels = names(cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])]))
    # ### the id column in the plot_df should be a factor
    # p[[3]] <- LabelClusters(plot = p[[3]], id = "direction2", col = "black")
    # 
    # ### save the plots
    # g <- arrangeGrob(grobs = p,
    #                  nrow = 3,
    #                  ncol = 1,
    #                  top = "")
    # ggsave(file = paste0(result_dir, paste(conds, collapse = "_vs_"),
    #                      "_RNAMagnet_Result_AD.png"), g, width = 20, height = 24, dpi = 300)
    
    ### draw a beeswarm plot with the adhesiveness info
    ggplot(plot_df, aes_string(x="cluster_color", y="adhesiveness")) +
      theme_classic(base_size = 16) +
      geom_boxplot() +
      geom_beeswarm(aes_string(color="direction"), na.rm = TRUE) +
      stat_compare_means() +
      labs(x = "", y = "Adhesiveness") +
      theme(legend.position="right")
    ggsave(file = paste0(result_dir, paste(conds, collapse = "_vs_"),
                         "_RNAMagnet_Beeswarm_AD.png"), width = 20, height = 12, dpi = 300)
    
    # ### draw a scatter plot with the specificity info
    # p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    #   geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
    #   xlab(x_lab) + ylab(x_lab) +
    #   labs(col="Cluster") +
    #   ggtitle(paste("UMAP with Cell Type")) +
    #   scale_color_brewer(palette="Dark2") +
    #   theme_classic(base_size = 16)
    # 
    # ### get colors for the clustering result
    # cell_colors_clust <- cell_pal(unique(SO@meta.data[,target_col]), hue_pal())
    # 
    # p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    #   geom_point(aes_string(col="direction", alpha="specificity"), size=2) +
    #   xlab(x_lab) + ylab(x_lab) +
    #   labs(col="Direction", alpha="Specificity Score") +
    #   ggtitle(paste("UMAP with Direction & Specificity Score")) +
    #   theme_classic(base_size = 16) +
    #   scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$direction)[order(unique(plot_df$direction))])],
    #                      labels = names(cell_colors_clust[as.character(unique(plot_df$direction)[order(unique(plot_df$direction))])]))
    # 
    # ### get colors for the clustering result
    # cell_colors_clust <- cell_pal(unique(SO@meta.data$Annotation), hue_pal())
    # 
    # plot_df$direction2 <- factor(plot_df$direction2, levels = unique(plot_df$direction2))
    # p[[3]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    #   geom_point(aes_string(col="direction2", alpha="specificity2"), size=2) +
    #   xlab(x_lab) + ylab(x_lab) +
    #   labs(col="Direction", alpha="Specificity") +
    #   ggtitle("UMAP with Direction & Specificity") +
    #   theme_classic(base_size = 16) +
    #   scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])],
    #                      labels = names(cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])]))
    # ### the id column in the plot_df should be a factor
    # p[[3]] <- LabelClusters(plot = p[[3]], id = "direction2", col = "black")
    # 
    # ### save the plots
    # g <- arrangeGrob(grobs = p,
    #                  nrow = 3,
    #                  ncol = 1,
    #                  top = "")
    # ggsave(file = paste0(result_dir, paste(conds, collapse = "_vs_"),
    #                      "_RNAMagnet_Result_SP.png"), g, width = 20, height = 24, dpi = 300)
    
    ### draw a beeswarm plot with the adhesiveness info
    ggplot(plot_df, aes_string(x="cluster_color", y="specificity")) +
      theme_classic(base_size = 16) +
      geom_boxplot() +
      geom_beeswarm(aes_string(color="direction"), na.rm = TRUE) +
      stat_compare_means() +
      labs(x = "", y = "Specificity") +
      theme(legend.position="right")
    ggsave(file = paste0(result_dir, paste(conds, collapse = "_vs_"),
                         "_RNAMagnet_Beeswarm_SP.png"), width = 20, height = 12, dpi = 300)
    
    ### run RNAMagnet signaling
    SO <- SetIdent(object = SO,
                   cells = rownames(SO@meta.data),
                   value = SO@meta.data[,target_col])
    result <- RNAMagnetSignaling(SO)
    
    ### get all the interaction list
    interaction_list <- NULL
    for(clust1 in conds) {
      for(clust2 in conds) {
        il <- getRNAMagnetGenes(result, clust1, clust2, thresh = 0)
        if(nrow(il) > 0) {
          il$ligand_cluster <- clust1
          il$receptor_cluster <- clust2
          if(is.null(interaction_list)) {
            interaction_list <- il
          } else {
            interaction_list <- rbind(interaction_list, il)
          }
        }
      }
    }
    interaction_list <- cbind(interaction_list[3:4], interaction_list[1:2])
    colnames(interaction_list) <- c("Ligand_Cluster", "Receptor_Cluster", "Interaction_Score", "Interaction_Pair")
    
    ### write out the interaction list
    write.xlsx2(interaction_list,
                file = paste0(result_dir, paste(conds, collapse = "_vs_"),
                              "_Signaling_RNAMagnet_Interaction_List.xlsx"),
                sheetName = paste0("Signaling_RNAMagnet_Interaction_List"),
                row.names = FALSE)
    
  }
  
  ### Adult LT-HSCs vs Adult CARs
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "Annotation2",
                               conds = c("LTHSC", "CARs"),
                               dim_method = "UMAP",
                               result_dir = outputDir2)
  
  ### Adult LT-HSCs vs Adult SECs
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "Annotation2",
                               conds = c("LTHSC", "SECs"),
                               dim_method = "UMAP",
                               result_dir = outputDir2)
  
  ### Adult LT-HSCs vs Adult AECs
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "Annotation2",
                               conds = c("LTHSC", "AECs"),
                               dim_method = "UMAP",
                               result_dir = outputDir2)
  
  ### Adult LT-HSCs vs Adult F1
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "Annotation2",
                               conds = c("LTHSC", "F-1"),
                               dim_method = "UMAP",
                               result_dir = outputDir2)
  
  ### Adult LT-HSCs vs Adult F4
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "Annotation2",
                               conds = c("LTHSC", "F-4"),
                               dim_method = "UMAP",
                               result_dir = outputDir2)
  
  ### Adult LT-HSCs vs Adult F5
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "Annotation2",
                               conds = c("LTHSC", "F-5"),
                               dim_method = "UMAP",
                               result_dir = outputDir2)
  
  
  ### RNAMagnet for the following:
  ### E16.5 LT-HSCs vs E16.5 Stroma
  ### E18.5 LT-HSCs vs E18.5 Stroma
  ### P0 LT-HSCs vs P0 Stroma
  
  ### set the ident of the object with the group info
  Updated_Seurat_Obj <- SetIdent(object = Updated_Seurat_Obj,
                                 cells = rownames(Updated_Seurat_Obj@meta.data),
                                 value = Updated_Seurat_Obj@meta.data$Development)
  
  ### E16.5 LT-HSCs vs E16.5 Stroma
  ### split the Seurat obj based on the given info
  Combined_Adult_Seurat_Obj <- subset(Updated_Seurat_Obj, idents="E16")
  
  ### check the UMAP for all the cells
  umap_plot <- DimPlot(Combined_Adult_Seurat_Obj, reduction = "umap", group.by = "Cell_Type", pt.size = 1.5) +
    labs(title = paste0("UMAP_E16"))
  umap_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  plot(umap_plot)
  
  ### run the RNAMagnet wrapper function
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "HSPC",
                               conds = c("LTHSC", "Stroma"),
                               dim_method = "PCA",
                               result_dir = paste0(outputDir2, "E16.5"))
  
  ### E18.5 LT-HSCs vs E18.5 Stroma
  ### split the Seurat obj based on the given info
  Combined_Adult_Seurat_Obj <- subset(Updated_Seurat_Obj, idents="E18")
  
  ### check the UMAP for all the cells
  umap_plot <- DimPlot(Combined_Adult_Seurat_Obj, reduction = "umap", group.by = "Cell_Type", pt.size = 1.5) +
    labs(title = paste0("UMAP_E16"))
  umap_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  plot(umap_plot)
  
  ### run the RNAMagnet wrapper function
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "HSPC",
                               conds = c("LTHSC", "Stroma"),
                               dim_method = "PCA",
                               result_dir = paste0(outputDir2, "E18.5"))
  
  ### P0 LT-HSCs vs P0 Stroma
  ### split the Seurat obj based on the given info
  Combined_Adult_Seurat_Obj <- subset(Updated_Seurat_Obj, idents="P0")
  
  ### check the UMAP for all the cells
  umap_plot <- DimPlot(Combined_Adult_Seurat_Obj, reduction = "umap", group.by = "Cell_Type", pt.size = 1.5) +
    labs(title = paste0("UMAP_E16"))
  umap_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  plot(umap_plot)
  
  ### run the RNAMagnet wrapper function
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "HSPC",
                               conds = c("LTHSC", "Stroma"),
                               dim_method = "PCA",
                               result_dir = paste0(outputDir2, "P0"))
  
  
  ### #3 analysis
  
  ### a function for the third analysis
  third_analysis <- function(Seurat_Object,
                             target_col,
                             comp1,
                             comp2,
                             time_point,
                             dim_method,
                             result_dir) {
    
    ### set output directory
    result_dir <- paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                         "_RNAMagnet_Result/", time_point, "/")
    dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
    
    ### set group info to the metadata
    Seurat_Object@meta.data$Group <- paste0(Seurat_Object@meta.data[,target_col], "_", Seurat_Object@meta.data$Development)
    
    ### all the comps
    comps <- union(paste0(comp1, "_", time_point), paste0(comp2, "_", time_point))
    
    ### set the ident of the object with the specified info
    Seurat_Object <- SetIdent(object = Seurat_Object,
                              cells = rownames(Seurat_Object@meta.data),
                              value = Seurat_Object@meta.data$Group)
    
    ### only keep the specified cells
    Seurat_Object <- subset(Seurat_Object, idents=comps)
    
    ### rownames in the meta.data should be in the same order as colnames in the counts
    Seurat_Object@meta.data <- Seurat_Object@meta.data[colnames(Seurat_Object@assays$RNA@counts),]
    
    ### preprocessing
    Seurat_Object <- FindVariableFeatures(Seurat_Object)
    Seurat_Object <- ScaleData(Seurat_Object)
    
    ### run PCA/UMAP
    if(dim_method == "PCA") {
      Seurat_Object <- RunPCA(Seurat_Object, npcs = 15)
      dim_map <- Embeddings(Seurat_Object, reduction = "pca")[rownames(Seurat_Object@meta.data),]
      x_lab <- "PC1"
      y_lab <- "PC2"
    } else if(dim_method == "UMAP") {
      Seurat_Object <- RunUMAP(Seurat_Object, dims = 1:15)
      dim_map <- Embeddings(Seurat_Object, reduction = "umap")[rownames(Seurat_Object@meta.data),]
      x_lab <- "UMAP1"
      y_lab <- "UMAP2"
    } else {
      stop("ERROR: dim_method not PCA nor UMAP.")
    }
    
    ### draw a UMAP with Trent's annotation
    Seurat_Object@meta.data$Annotation <- factor(Seurat_Object@meta.data$Annotation)
    umap_plot <- DimPlot(Seurat_Object, reduction = "umap", group.by = "Annotation", pt.size = 1.5) +
      labs(title = paste0("UMAP_With_the_Annotation"))
    umap_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.5
    umap_plot <- LabelClusters(plot = umap_plot, id = "Annotation", col = "black")
    ggsave(file = paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                         "_UMAP_with_the_annotation_", time_point, ".png"),
           width = 15, height = 10, dpi = 300)
    
    ### run RNAMagnet with anchors
    ### Warning: the Idents(Seurat_Object) should be along with the given 'anchors' input
    ### Another ERROR:
    ### The following code is fixed due to an ERROR: 'drop=FALSE'
    #compute mean gene expression level per population
    # out@anchors <- do.call(cbind, lapply(anchors, function(id) {
    #   apply(resolvedRawData[,seurat.ident == id,drop=FALSE],1,mean)
    # }));
    # colnames(out@anchors) <- anchors
    Seurat_Object <- SetIdent(object = Seurat_Object,
                              cells = rownames(Seurat_Object@meta.data),
                              value = Seurat_Object@meta.data$Group)
    result <- RNAMagnetAnchors(Seurat_Object,
                               anchors = unique(Seurat_Object@meta.data$Group))
    Seurat_Object <- SetIdent(object = Seurat_Object,
                              cells = rownames(Seurat_Object@meta.data),
                              value = Seurat_Object@meta.data$Annotation)
    result2 <- RNAMagnetAnchors(Seurat_Object,
                                anchors = unique(Seurat_Object@meta.data$Annotation))
    
    ### add RNAMagnet info to the seurat object
    Seurat_Object@meta.data$direction <- as.character(result[rownames(Seurat_Object@meta.data),"direction"])
    Seurat_Object@meta.data$direction2 <- as.character(result2[rownames(Seurat_Object@meta.data),"direction"])
    Seurat_Object@meta.data$adhesiveness <- as.numeric(result[rownames(Seurat_Object@meta.data),"adhesiveness"])
    Seurat_Object@meta.data$adhesiveness2 <- as.numeric(result2[rownames(Seurat_Object@meta.data),"adhesiveness"])
    Seurat_Object@meta.data$specificity <- as.numeric(sapply(rownames(Seurat_Object@meta.data), function(x) {
      result[x, result[x,"direction"]]
    }))
    Seurat_Object@meta.data$specificity2 <- as.numeric(sapply(rownames(Seurat_Object@meta.data), function(x) {
      result2[x, result2[x,"direction"]]
    }))
    
    ### make a data frame for ggplot
    plot_df <- data.frame(X=dim_map[rownames(Seurat_Object@meta.data),1],
                          Y=dim_map[rownames(Seurat_Object@meta.data),2],
                          cluster_color = Seurat_Object@meta.data$Group,
                          annotation_color = Seurat_Object@meta.data$Annotation,
                          direction = Seurat_Object@meta.data$direction,
                          direction2 = Seurat_Object@meta.data$direction2,
                          adhesiveness = Seurat_Object@meta.data$adhesiveness,
                          adhesiveness2 = Seurat_Object@meta.data$adhesiveness2,
                          specificity = Seurat_Object@meta.data$specificity,
                          specificity2 = Seurat_Object@meta.data$specificity2,
                          stringsAsFactors = FALSE, check.names = FALSE)
    
    ### get colors for the clustering result
    cell_colors_clust <- cell_pal(unique(Seurat_Object@meta.data$Annotation), hue_pal())
    
    ### scatter plot
    p <- list()
    
    ### draw a scatter plot with the adhesiveness info
    p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
      xlab(x_lab) + ylab(y_lab) +
      labs(col="Cluster") +
      ggtitle("UMAP with Cell Type") +
      scale_color_brewer(palette="Dark2") +
      theme_classic(base_size = 16)
    
    p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="annotation_color"), size=2, alpha=0.5) +
      xlab(x_lab) + ylab(y_lab) +
      labs(col="Cluster") +
      ggtitle("UMAP with Cell Type") +
      theme_classic(base_size = 16)
    ### the id column in the plot_df should be a factor
    p[[2]] <- LabelClusters(plot = p[[2]], id = "annotation_color", col = "black")
    
    p[[3]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="direction", alpha="adhesiveness"), size=2) +
      xlab(x_lab) + ylab(y_lab) +
      labs(col="Direction", alpha="Adhesiveness") +
      ggtitle("UMAP with Direction & Adhesiveness") +
      scale_color_brewer(palette="Set1") +
      theme_classic(base_size = 16)
    
    plot_df$direction2 <- factor(plot_df$direction2, levels = unique(plot_df$direction2))
    p[[4]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="direction2", alpha="adhesiveness2"), size=2) +
      xlab(x_lab) + ylab(y_lab) +
      labs(col="Direction", alpha="Adhesiveness") +
      ggtitle("UMAP with Direction & Adhesiveness") +
      theme_classic(base_size = 16) +
      scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])],
                         labels = names(cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])]))
    ### the id column in the plot_df should be a factor
    p[[4]] <- LabelClusters(plot = p[[4]], id = "direction2", col = "black")
    
    ### save the plots
    g <- arrangeGrob(grobs = p,
                     nrow = 2,
                     ncol = 2,
                     top = "")
    ggsave(file = paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                         "_RNAMagnet_Result_AD_", time_point, "_New_Annotation.png"), g, width = 28, height = 20, dpi = 300)
    
    ### draw a scatter plot with the specificity info
    p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
      xlab(x_lab) + ylab(y_lab) +
      labs(col="Cluster") +
      ggtitle("UMAP with Cell Type") +
      scale_color_brewer(palette="Dark2") +
      theme_classic(base_size = 16)
    
    p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="annotation_color"), size=2, alpha=0.5) +
      xlab(x_lab) + ylab(y_lab) +
      labs(col="Cluster") +
      ggtitle("UMAP with Cell Type") +
      theme_classic(base_size = 16)
    ### the id column in the plot_df should be a factor
    p[[2]] <- LabelClusters(plot = p[[2]], id = "annotation_color", col = "black")
    
    p[[3]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="direction", alpha="specificity"), size=2) +
      xlab(x_lab) + ylab(y_lab) +
      labs(col="Direction", alpha="Specificity") +
      ggtitle("UMAP with Direction & Specificity") +
      scale_color_brewer(palette="Set1") +
      theme_classic(base_size = 16)
    
    plot_df$direction2 <- factor(plot_df$direction2, levels = unique(plot_df$direction2))
    p[[4]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="direction2", alpha="specificity2"), size=2) +
      xlab(x_lab) + ylab(y_lab) +
      labs(col="Direction", alpha="Specificity") +
      ggtitle("UMAP with Direction & Specificity") +
      theme_classic(base_size = 16) +
      scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])],
                         labels = names(cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])]))
    ### the id column in the plot_df should be a factor
    p[[4]] <- LabelClusters(plot = p[[4]], id = "direction2", col = "black")
    
    ### save the plots
    g <- arrangeGrob(grobs = p,
                     nrow = 2,
                     ncol = 2,
                     top = "")
    ggsave(file = paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                         "_RNAMagnet_Result_SP_", time_point, "_New_Annotation.png"), g, width = 28, height = 20, dpi = 300)
    
  }
  
  ### adult MPP2/3/4/ST-HSCs vs adult stroma
  third_analysis(Seurat_Object = Updated_Seurat_Obj,
                 target_col = "HSPC",
                 comp1 = c("STHSC", "MPP2", "MPP3", "MPP4"),
                 comp2 = "Stroma",
                 time_point = "ADULT",
                 dim_method = "UMAP",
                 result_dir=outputDir2)
  
  ### E16.5 MPP2/3/4/ST-HSCs vs E16.5 Stroma
  third_analysis(Seurat_Object = Updated_Seurat_Obj,
                 target_col = "HSPC",
                 comp1 = c("STHSC", "MPP2", "MPP3", "MPP4"),
                 comp2 = "Stroma",
                 time_point = "E16",
                 dim_method = "UMAP",
                 result_dir=outputDir2)
  
  ### E18.5 MPP2/3/4/ST-HSCs vs E18.5 Stroma
  third_analysis(Seurat_Object = Updated_Seurat_Obj,
                 target_col = "HSPC",
                 comp1 = c("STHSC", "MPP2", "MPP3", "MPP4"),
                 comp2 = "Stroma",
                 time_point = "E18",
                 dim_method = "UMAP",
                 result_dir=outputDir2)
  
  ### P0 MPP2/3/4/ST-HSCs vs P0 stroma
  third_analysis(Seurat_Object = Updated_Seurat_Obj,
                 target_col = "HSPC",
                 comp1 = c("STHSC", "MPP2", "MPP3", "MPP4"),
                 comp2 = "Stroma",
                 time_point = "P0",
                 dim_method = "UMAP",
                 result_dir=outputDir2)
  
}

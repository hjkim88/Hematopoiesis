###
#   File name : RNAMagnet_Analysis.R
#   Author    : Hyunjin Kim
#   Date      : Aug 24, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : 1. For Heme objects E16, E18, P0, and Adult, run RNA-Magnet for physical interactions of
#                  each cluster within the object and against the stroma object of the same timepoint.
#                  Also, run RNA-Magnet for physical interactions of each cluster within each stroma object.
#                  * What Trent would like to see is:
#                   a) A tSNE-plot showing the "direction", "adhesiveness", and Specificity Score of the cells
#                      in each Heme object for all anchor populations (Heme and Stroma) at the same timepoint.
#                   b) The same data in a heat map.
#               2. For Heme objects E16, E18, P0, and Adult, run RNA-Magnet to infer the putative signaling
#                  interactions between each cluster within the object and between the stroma objects of the same
#                  timepoint. Also, run RNA-Magnet for putative signaling interactions bewteen cluster within
#                  each stroma object.
#                  * What Trent would like to see here is:
#                   a) A summary of the "specificity" score at the cluster level i) within each Heme object,
#                      ii) between the Heme object and the Stroma object at the same timepoint, and
#                      iii) within each Stroma object, depicted as a signaling network.
#                   b) A list of the molecules within each cluster mediating these interactions.
#                   c) A tSNE-Plot showing the "specificity" of the cells in each Heme object for
#                      i) each cluster within the object and ii) each cluster within the stroma object of the
#                      same timepoint, as well as a tSNE-Plot showing this same value in each Stroma object
#                      for each cluster within the object.
#   
#   Instruction
#               1. Source("RNAMagnet_Analysis.R")
#               2. Run the function "rna_magnet_trent" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_RNAMagnet_Analysis.R/RNAMagnet_Analysis.R")
#               > rna_magnet_trent(Robj_path="C:/Users/hkim8/SJ/HSPC Analyses/Stroma and Heme/",
#                                  Seurat_RObj_path="./data/Combined_BM_Seurat_Obj.RDATA",
#                                  outputDir="./results/RNA-Magnet/")
###

rna_magnet_trent <- function(Robj_path="C:/Users/hkim8/SJ/HSPC Analyses/Stroma and Heme/",
                             Seurat_RObj_path="./data/Combined_BM_Seurat_Obj.RDATA",
                             outputDir="./results/RNA-Magnet/") {
  
  ### load libraries
  options(java.parameters = "-Xmx10240m")
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(remotes, quietly = TRUE)) {
    install.packages("remotes")
    require(remotes, quietly = TRUE)
  }
  if(!require(RNAMagnet, quietly = TRUE)) {
    remotes::install_github("veltenlab/rnamagnet")
    require(RNAMagnet, quietly = TRUE)
  }
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(reticulate, quietly = TRUE)) {
    install.packages("reticulate")
    require(reticulate, quietly = TRUE)
  }
  if(!require(Rmagic, quietly = TRUE)) {
    install.packages("Rmagic")
    require(Rmagic, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(scales, quietly = TRUE)) {
    install.packages("scales")
    library(scales, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  
  ### see python environment since RNAMagnet uses the python module 'magic'
  # conda_create("r-reticulate")
  # conda_install("r-reticulate", "python-magic")
  use_condaenv("r-reticulate")
  py_config()
  py_module_available("magic")
  
  ### get Robj file list
  f_list <- list.files(Robj_path, pattern = ".Robj$")
  
  ### determine necessary variables
  time_points <- c("E16", "E18", "P0", "ADULT")
  cell_type <- c("Heme", "Stroma")
  possible_names <- as.vector(sapply(cell_type, function(x) paste0(x, time_points)))
  possible_names_mat <- data.frame(Names=possible_names,
                                   Type=as.vector(sapply(cell_type, function(x) rep(x, length(time_points)))),
                                   Time=rep(time_points, length(cell_type)),
                                   stringsAsFactors = FALSE, check.names = FALSE)
  row.names(possible_names_mat) <- possible_names
  
  ### keep the required files only
  f_list <- f_list[which(f_list %in% paste0(possible_names, "_regress.Robj"))]
  
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
  
  ### for each of the Robj file, run the RNAMagnet
  Combined_Seurat_Obj <- NULL
  for(f in f_list) {
    
    ### load the Seurat object and save the object name
    tmp_env <- new.env()
    load(paste0(Robj_path, f), tmp_env)
    obj_name <- ls(tmp_env)
    assign("Seurat_Obj", get(obj_name, envir = tmp_env))
    rm(tmp_env)
    gc()
    
    ### rownames in the meta.data should be in the same order as colnames in the counts
    Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
    
    ### active assay = "RNA"
    Seurat_Obj@active.assay <- "RNA"
    
    ### set output directory
    outputDir2 <- paste0(outputDir, f, "/")
    dir.create(outputDir2, recursive = TRUE, showWarnings = FALSE)
    
    ### annotate the file name to the R object
    tissue_name <- strsplit(f, split = "_", fixed = TRUE)[[1]][1]
    Seurat_Obj@meta.data$Tissue <- tissue_name
    
    ### annotate the development
    Seurat_Obj@meta.data$Development <- possible_names_mat[tissue_name, "Time"]
    
    ### annotate the Cell type
    Seurat_Obj@meta.data$Cell_Type <- possible_names_mat[tissue_name, "Type"]
    
    ### combine the Seurat R objects
    if(is.null(Combined_Seurat_Obj)) {
      Combined_Seurat_Obj <- Seurat_Obj
    } else {
      Combined_Seurat_Obj <- merge(x = Combined_Seurat_Obj, y = Seurat_Obj,
                                   project = Seurat_Obj@project.name)
    }
    
    ### run only if there are enough number of cells in the R object
    if(nrow(Seurat_Obj@meta.data) > 1 && identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data))) {
      
      ### clustering on the R object
      Seurat_Obj <- FindNeighbors(Seurat_Obj, dims = 1:5, k.param = 5)
      Seurat_Obj <- FindClusters(Seurat_Obj, resolution = 0.4)
      Seurat_Obj@meta.data$new_clusts <- Idents(Seurat_Obj)
      
      #
      ### run RNAMagnet
      #
      ### get PCA matrix
      pca_map <- Embeddings(Seurat_Obj, reduction = "pca")[rownames(Seurat_Obj@meta.data),1:5]
      
      ### set the ident of the seurat object with the cluster info
      Seurat_Obj <- SetIdent(object = Seurat_Obj,
                             cells = rownames(Seurat_Obj@meta.data),
                             value = Seurat_Obj@meta.data$new_clusts)
      
      ### run RNAMagnet
      result <- RNAMagnetAnchors(Seurat_Obj,
                                 anchors = levels(Seurat_Obj@meta.data$new_clusts))
      
      ### make file name
      obj <- strsplit(f, split = "_", fixed = TRUE)[[1]][1]
      
      ### write the result as an Excel file
      write.xlsx2(data.frame(Cell=rownames(result), result,
                             stringsAsFactors = FALSE, check.names = FALSE),
                  file = paste0(outputDir2, "RNAMagnet_Result_", obj, ".xlsx"),
                  sheetName = paste0("RNAMagnet_", obj),
                  row.names = FALSE)
      
      ### add RNAMagnet info to the seurat object
      Seurat_Obj@meta.data$direction <- as.character(result[rownames(Seurat_Obj@meta.data),"direction"])
      Seurat_Obj@meta.data$adhesiveness <- as.numeric(result[rownames(Seurat_Obj@meta.data),"adhesiveness"])
      Seurat_Obj@meta.data$specificity <- as.numeric(sapply(rownames(Seurat_Obj@meta.data), function(x) {
        result[x, paste0("X", result[x,"direction"])]
      }))
      
      ### check the order
      print(identical(rownames(Seurat_Obj@meta.data), rownames(pca_map)))
      
      ### make a data frame for ggplot
      plot_df <- data.frame(X=pca_map[rownames(Seurat_Obj@meta.data),"PC_1"],
                            Y=pca_map[rownames(Seurat_Obj@meta.data),"PC_2"],
                            group_color = Seurat_Obj@meta.data$direction,
                            group_alpha = Seurat_Obj@meta.data$adhesiveness,
                            cluster_color = Seurat_Obj@meta.data$new_clusts,
                            specificity = Seurat_Obj@meta.data$specificity,
                            stringsAsFactors = FALSE, check.names = FALSE)
      
      ### factorize the direction columns
      plot_df$group_color <- factor(plot_df$group_color, levels = levels(plot_df$cluster_color)) 
      
      ### get colors for the clustering result
      cell_colors_clust <- cell_pal(levels(Seurat_Obj@meta.data$new_clusts), hue_pal())
      
      ### scatter plot
      p <- list()
      
      ### original PCA
      p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
        geom_point(aes_string(col="cluster_color"), size=2) +
        xlab("PC1") + ylab("PC2") +
        labs(col="Cluster") +
        ggtitle("PCA with Cluster Info") +
        theme_classic(base_size = 16) +
        scale_color_manual(values = cell_colors_clust, labels = names(cell_colors_clust))
      p[[1]] <- LabelClusters(plot = p[[1]], id = "cluster_color")
      
      ### draw a scatter plot with the adhesiveness info
      p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
        geom_point(aes_string(col="group_color", alpha="group_alpha"), size=2) +
        xlab("PC1") + ylab("PC2") +
        labs(col="Direction", alpha="Adhesiveness") +
        ggtitle("PCA with Direction & Adhesiveness") +
        theme_classic(base_size = 16) +
        scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$group_color)[order(unique(plot_df$group_color))])],
                           labels = names(cell_colors_clust[as.character(unique(plot_df$group_color)[order(unique(plot_df$group_color))])]))
      p[[2]] <- LabelClusters(plot = p[[2]], id = "cluster_color")
      
      ### save the plots
      g <- arrangeGrob(grobs = p,
                       nrow = 2,
                       ncol = 1,
                       top = "")
      ggsave(file = paste0(outputDir2, "PCA_RNAMagnet_Result_AD_", obj, ".png"), g, width = 20, height = 12, dpi = 300)
      
      ### draw a scatter plot with the specificity info
      p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
        geom_point(aes_string(col="group_color", alpha="specificity"), size=2) +
        xlab("PC1") + ylab("PC2") +
        labs(col="Direction", alpha="Specificity Score") +
        ggtitle("PCA with Direction & Specificity Score") +
        theme_classic(base_size = 16) +
        scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$group_color)[order(unique(plot_df$group_color))])],
                           labels = names(cell_colors_clust[as.character(unique(plot_df$group_color)[order(unique(plot_df$group_color))])]))
      p[[2]] <- LabelClusters(plot = p[[2]], id = "cluster_color")
      
      ### save the plots
      g <- arrangeGrob(grobs = p,
                       nrow = 2,
                       ncol = 1,
                       top = "")
      ggsave(file = paste0(outputDir2, "PCA_RNAMagnet_Result_SP_", obj, ".png"), g, width = 20, height = 12, dpi = 300)
      
      #
      ### heatmap with specificity scores - col side color bars with direction & adhesiveness
      #
      
      
      
      
    }
    
  }
  
}

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
  
  ### result directory
  outputDir2 <- paste0(outputDir, "RNA_Magnet/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### see python environment since RNAMagnet uses the python module 'magic'
  # conda_create("r-reticulate")
  # conda_install("r-reticulate", "python-magic")
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
  
  ### simply run RNAMagnet only to get signaling interactions
  ### SO: Seurat object
  ### target_col: Target column name of the given seurat object's meta.data
  ### conds: a character vector of all the conditions in the 'target_col' that will be used
  ### result_dir: output directory
  simple_RNAMagnet_interaction <- function(SO, target_col, conds, result_dir) {
    
    result_dir <- paste0(result_dir, "/", paste(conds, collapse = "_vs_"), "/")
    dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
    
    ### set the ident of the object with the HSPC type
    SO <- SetIdent(object = SO,
                   cells = rownames(SO@meta.data),
                   value = SO@meta.data[,target_col])
    
    ### split the Seurat obj based on the given info
    SO <- subset(SO, idents=conds)
    
    ### run UMAP
    SO <- RunUMAP(SO, dims = 1:15)
    dim_map <- Embeddings(SO, reduction = "umap")[rownames(SO@meta.data),]
    
    ### run RNAMagnet with anchors
    result <- RNAMagnetAnchors(SO,
                               anchors = unique(SO@meta.data[,target_col]))
    
    ### write the result as an Excel file
    write.xlsx2(data.frame(Cell=rownames(result), result,
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(result_dir, paste(conds, collapse = "_vs_"),
                              "_RNAMagnet_Result.xlsx"),
                sheetName = paste0("RNAMagnet_Result"),
                row.names = FALSE)
    
    ### add RNAMagnet info to the seurat object
    SO@meta.data$direction <- as.character(result[rownames(SO@meta.data),"direction"])
    SO@meta.data$adhesiveness <- as.numeric(result[rownames(SO@meta.data),"adhesiveness"])
    SO@meta.data$specificity <- as.numeric(sapply(rownames(SO@meta.data), function(x) {
      result[x, result[x,"direction"]]
    }))
    
    ### make a data frame for ggplot
    plot_df <- data.frame(X=dim_map[rownames(SO@meta.data),1],
                          Y=dim_map[rownames(SO@meta.data),2],
                          direction = SO@meta.data$direction,
                          adhesiveness = SO@meta.data$adhesiveness,
                          cluster_color = SO@meta.data[,target_col],
                          specificity = SO@meta.data$specificity,
                          stringsAsFactors = FALSE, check.names = FALSE)
    
    ### get colors for the clustering result
    cell_colors_clust <- cell_pal(unique(SO@meta.data[,target_col]), hue_pal())
    
    ### scatter plot
    p <- list()
    
    ### draw a scatter plot with the adhesiveness info
    p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Cluster") +
      ggtitle(paste0("UMAP with Cell Type")) +
      scale_color_brewer(palette="Dark2") +
      theme_classic(base_size = 16)
    
    p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="direction", alpha="adhesiveness"), size=2) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Direction", alpha="Adhesiveness") +
      ggtitle(paste("UMAP with Direction & Adhesiveness")) +
      theme_classic(base_size = 16) +
      scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$direction)[order(unique(plot_df$direction))])],
                         labels = names(cell_colors_clust[as.character(unique(plot_df$direction)[order(unique(plot_df$direction))])]))
    
    ### save the plots
    g <- arrangeGrob(grobs = p,
                     nrow = 2,
                     ncol = 1,
                     top = "")
    ggsave(file = paste0(result_dir, paste(conds, collapse = "_vs_"),
                         "_RNAMagnet_Result_AD.png"), g, width = 20, height = 12, dpi = 300)
    
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
    
    ### draw a scatter plot with the specificity info
    p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Cluster") +
      ggtitle(paste("UMAP with Cell Type")) +
      scale_color_brewer(palette="Dark2") +
      theme_classic(base_size = 16)
    
    p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="direction", alpha="specificity"), size=2) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Direction", alpha="Specificity Score") +
      ggtitle(paste("UMAP with Direction & Specificity Score")) +
      theme_classic(base_size = 16) +
      scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$direction)[order(unique(plot_df$direction))])],
                         labels = names(cell_colors_clust[as.character(unique(plot_df$direction)[order(unique(plot_df$direction))])]))
    
    ### save the plots
    g <- arrangeGrob(grobs = p,
                     nrow = 2,
                     ncol = 1,
                     top = "")
    ggsave(file = paste0(result_dir, paste(conds, collapse = "_vs_"),
                         "_RNAMagnet_Result_SP.png"), g, width = 20, height = 12, dpi = 300)
    
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
                               result_dir = outputDir2)
  
  ### Adult LT-HSCs vs Adult SECs
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "Annotation2",
                               conds = c("LTHSC", "SECs"),
                               result_dir = outputDir2)
  
  ### Adult LT-HSCs vs Adult AECs
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "Annotation2",
                               conds = c("LTHSC", "AECs"),
                               result_dir = outputDir2)
  
  ### Adult LT-HSCs vs Adult F1
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "Annotation2",
                               conds = c("LTHSC", "F-1"),
                               result_dir = outputDir2)
  
  ### Adult LT-HSCs vs Adult F4
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "Annotation2",
                               conds = c("LTHSC", "F-4"),
                               result_dir = outputDir2)
  
  ### Adult LT-HSCs vs Adult F5
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "Annotation2",
                               conds = c("LTHSC", "F-5"),
                               result_dir = outputDir2)
  
  
  ### RNAMagnet for the following:
  ### E16.5 LT-HSCs vs E16.5 Stroma
  ### E18.5 LT-HSCs vs E18.5 Stroma
  ### P0 LT-HSCs vs P0 Stroma
  
  ### E16.5 LT-HSCs vs E16.5 Stroma
  ### split the Seurat obj based on the given info
  Combined_Adult_Seurat_Obj <- subset(Updated_Seurat_Obj, idents="E16")
  
  ### run the RNAMagnet wrapper function
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "HSPC",
                               conds = c("LTHSC", "Stroma"),
                               result_dir = paste0(outputDir2, "E16.5"))
  
  ### E18.5 LT-HSCs vs E18.5 Stroma
  ### split the Seurat obj based on the given info
  Combined_Adult_Seurat_Obj <- subset(Updated_Seurat_Obj, idents="E18")
  
  ### run the RNAMagnet wrapper function
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "HSPC",
                               conds = c("LTHSC", "Stroma"),
                               result_dir = paste0(outputDir2, "E18.5"))
  
  ### P0 LT-HSCs vs P0 Stroma
  ### split the Seurat obj based on the given info
  Combined_Adult_Seurat_Obj <- subset(Updated_Seurat_Obj, idents="P0")
  
  ### run the RNAMagnet wrapper function
  simple_RNAMagnet_interaction(SO = Combined_Adult_Seurat_Obj,
                               target_col = "HSPC",
                               conds = c("LTHSC", "Stroma"),
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
    Seurat_Object@meta.data$Group <- paste0(Seurat_Object@meta.data[,target_col], "_", Seurat_Object@meta.data$Time)
    
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
    } else if(dim_method == "UMAP") {
      Seurat_Object <- RunUMAP(Seurat_Object, dims = 1:15)
      dim_map <- Embeddings(Seurat_Object, reduction = "umap")[rownames(Seurat_Object@meta.data),]
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
                          direction = Seurat_Object@meta.data$direction,
                          direction2 = Seurat_Object@meta.data$direction2,
                          adhesiveness = Seurat_Object@meta.data$adhesiveness,
                          adhesiveness2 = Seurat_Object@meta.data$adhesiveness2,
                          cluster_color = Seurat_Object@meta.data$Group,
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
      xlab("PC1") + ylab("PC2") +
      labs(col="Cluster") +
      ggtitle("UMAP with Cell Type") +
      scale_color_brewer(palette="Dark2") +
      theme_classic(base_size = 16)
    
    p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="direction", alpha="adhesiveness"), size=2) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Direction", alpha="Adhesiveness") +
      ggtitle("UMAP with Direction & Adhesiveness") +
      scale_color_brewer(palette="Set1") +
      theme_classic(base_size = 16)
    
    plot_df$direction2 <- factor(plot_df$direction2, levels = unique(plot_df$direction2))
    p[[3]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="direction2", alpha="adhesiveness2"), size=2) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Direction", alpha="Adhesiveness") +
      ggtitle("UMAP with Direction & Adhesiveness") +
      theme_classic(base_size = 16) +
      scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])],
                         labels = names(cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])]))
    ### the id column in the plot_df should be a factor
    p[[3]] <- LabelClusters(plot = p[[3]], id = "direction2", col = "black")
    
    ### save the plots
    g <- arrangeGrob(grobs = p,
                     nrow = 3,
                     ncol = 1,
                     top = "")
    ggsave(file = paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                         "_RNAMagnet_Result_AD_", time_point, "_New_Annotation.png"), g, width = 20, height = 20, dpi = 300)
    
    ### draw a scatter plot with the specificity info
    p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Cluster") +
      ggtitle("UMAP with Cell Type") +
      scale_color_brewer(palette="Dark2") +
      theme_classic(base_size = 16)
    
    p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="direction", alpha="specificity"), size=2) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Direction", alpha="Specificity") +
      ggtitle("UMAP with Direction & Specificity") +
      scale_color_brewer(palette="Set1") +
      theme_classic(base_size = 16)
    
    plot_df$direction2 <- factor(plot_df$direction2, levels = unique(plot_df$direction2))
    p[[3]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="direction2", alpha="specificity2"), size=2) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Direction", alpha="Specificity") +
      ggtitle("UMAP with Direction & Specificity") +
      theme_classic(base_size = 16) +
      scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])],
                         labels = names(cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])]))
    ### the id column in the plot_df should be a factor
    p[[3]] <- LabelClusters(plot = p[[3]], id = "direction2", col = "black")
    
    ### save the plots
    g <- arrangeGrob(grobs = p,
                     nrow = 3,
                     ncol = 1,
                     top = "")
    ggsave(file = paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                         "_RNAMagnet_Result_SP_", time_point, "_New_Annotation.png"), g, width = 20, height = 20, dpi = 300)
    
  }
  
  ### adult MPP2/3/4/ST-HSCs vs adult stroma
  third_analysis(Seurat_Object = Updated_Seurat_Obj,
                 target_col = "HSPC",
                 comp1 = c("STHSC", "MPP2", "MPP3", "MPP4"),
                 comp2 = "Stroma",
                 time_point = "Adult",
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

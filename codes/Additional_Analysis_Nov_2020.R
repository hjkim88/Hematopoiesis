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
  
  ### simply run RNAMagnet only to get signaling interactions
  ### SO: Seurat object
  ### target_col: Target column name of the given seurat object's meta.data
  ### conds: a character vector of all the conditions in the 'target_col' that will be used
  ### result_dir: output directory
  simple_RNAMagnet_interaction <- function(SO, target_col, conds, result_dir) {
    
    ### set the ident of the object with the HSPC type
    SO <- SetIdent(object = SO,
                   cells = rownames(SO@meta.data),
                   value = SO@meta.data[,target_col])
    
    ### split the Seurat obj based on the given info
    SO <- subset(SO, idents=conds)
    
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
    result_dir <- paste0(result_dir, "/", paste(conds, collapse = "_vs_"), "/")
    dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
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
  
}

###
#   File name : Additional_Analysis_Nov_2020.R
#   Author    : Hyunjin Kim
#   Date      : Nov 18, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Perform additional analyses for Trent's project based on the
#               "20201117 Meeting.txt".
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
  Updated_Seurat_Obj <- readRDS(file = Robj2_path)
  
  ### create new column for the analysis
  Updated_Seurat_Obj@meta.data$Dev_Anno <- paste0(Updated_Seurat_Obj@meta.data$Development, "_",
                                                  Updated_Seurat_Obj@meta.data$Annotation)
  
  ### determine necessary variables
  time_points <- c("E16", "E18", "P0", "ADULT")
  HSPC_populations <- c("LTHSC", "STHSC", "MPP2", "MPP3", "MPP4")
  
  ### set the ident of the object with the HSPC type
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
  ggsave(file = paste0(outputDir2, "UMAP_Combined_Adult_With_the_Annotation.png"), width = 15, height = 10, dpi = 300)
  
  ### Pull out stroma cells that are direct to LTHSCs/Other HSPCs and run RNA-Magnet - do the interaction scores change?
  
  
  
  
  
  
}

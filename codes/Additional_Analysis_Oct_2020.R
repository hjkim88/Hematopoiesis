###
#   File name : Additional_Analysis_Oct_2020.R
#   Author    : Hyunjin Kim
#   Date      : Oct 28, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Perform additional analyses for Trent's project based on the
#               "20201012 Analyses for Hyunjin Overview.xlsx".
#   
#   Instruction
#               1. Source("Additional_Analysis_Oct_2020.R")
#               2. Run the function "additional_analysis" - specify the input paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Additional_Analysis_Oct_2020.R/Additional_Analysis_Oct_2020.R")
#               > additional_analysis(Robj1_path="./data/Combined_Seurat_Obj.RDATA",
#                                     Robj2_path="./data/Combined_BM_Seurat_Obj.RDATA",
#                                     outputDir="./results/Additional_Oct2020/")
###

additional_analysis <- function(Robj1_path="./data/Combined_Seurat_Obj.RDATA",
                                Robj2_path="./data/Combined_BM_Seurat_Obj.RDATA",
                                outputDir="./results/Additional_Oct2020/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(slingshot, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("slingshot")
    require(slingshot, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  
  ### load the Seurat object and save the object name
  tmp_env <- new.env()
  load(paste0(Robj_path), tmp_env)
  obj_name <- ls(tmp_env)
  assign("Combined_Seurat_Obj", get(obj_name, envir = tmp_env))
  rm(tmp_env)
  gc()
  
  ### run PCA
  Combined_Seurat_Obj <- FindVariableFeatures(Combined_Seurat_Obj)
  Combined_Seurat_Obj <- ScaleData(Combined_Seurat_Obj)
  Combined_Seurat_Obj <- RunPCA(Combined_Seurat_Obj, npcs = 10)
  
  ### draw a PCA
  DimPlot(Combined_Seurat_Obj, reduction = "pca", group.by = "Tissue", pt.size = 1.5) +
    labs(title = paste0("PCA_Combined_Tissue"))
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Combined_Seurat_Obj)), rownames(Combined_Seurat_Obj@meta.data)))
  
  ### determine necessary variables
  time_points <- c("E16", "E18", "P0", "ADULT")
  HSPC_populations <- c("LTHSC", "STHSC", "MPP2", "MPP3", "MPP4")
  possible_names <- as.vector(sapply(time_points, function(x) paste0(x, HSPC_populations)))
  possible_names_mat <- data.frame(Names=possible_names,
                                   Time=as.vector(sapply(time_points, function(x) rep(x, length(HSPC_populations)))),
                                   HSPC=rep(HSPC_populations, length(time_points)),
                                   stringsAsFactors = FALSE, check.names = FALSE)
  row.names(possible_names_mat) <- possible_names
  
  #
  ### separate the combined into different cell groups
  #
  
  ### set the ident of the object with the HSPC type
  Combined_Seurat_Obj <- SetIdent(object = Combined_Seurat_Obj,
                                  cells = rownames(Combined_Seurat_Obj@meta.data),
                                  value = Combined_Seurat_Obj@meta.data$HSPC)
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Combined_Seurat_Obj)), rownames(Combined_Seurat_Obj@meta.data)))
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Combined_Seurat_Obj@meta.data <- Combined_Seurat_Obj@meta.data[colnames(Combined_Seurat_Obj@assays$RNA@counts),]
  
  ### active assay = "RNA"
  Combined_Seurat_Obj@active.assay <- "RNA"
  
  
  
  
  
}

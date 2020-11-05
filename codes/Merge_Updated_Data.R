###
#   File name : Merge_Updated_Data.R
#   Author    : Hyunjin Kim
#   Date      : Nov 4, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Combine the newly updated data and make a RDS file for future convenience
#   
#   Instruction
#               1. Source("Merge_Updated_Data.R")
#               2. Run the function "merge_data" - specify the input directory and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Merge_Updated_Data.R/Merge_Updated_Data.R")
#               > merge_data(Robj_path="C:/Users/hkim8/SJ/HSPC Analyses/Stroma and Heme/",
#                            Robj2_path="C:/Users/hkim8/SJ/HSPC Analyses/HSPC Subsets/",
#                            outputDir="./data/")
###

merge_data <- function(Robj_path="C:/Users/hkim8/SJ/HSPC Analyses/Stroma and Heme/",
                       Robj2_path="C:/Users/hkim8/SJ/HSPC Analyses/HSPC Subsets/",
                       outputDir="./data/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  
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
  f_list <- f_list[which(f_list %in% paste0(possible_names, "_regress2.Robj"))]
  
  ### for each of the Robj file, run the pseudotime analysis and combine the R objects for later use
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
    
  }
  
  
  ### get HSPC-specific Robj file list
  f_list <- list.files(Robj2_path, pattern = ".Robj$")
  
  ### determine necessary variables
  time_points <- c("E16", "E18", "P0", "ADULT")
  HSPC_populations <- c("LTHSC", "STHSC", "MPP2", "MPP3", "MPP4")
  possible_names <- as.vector(sapply(time_points, function(x) paste0(x, HSPC_populations)))
  possible_names_mat <- data.frame(Names=possible_names,
                                   Time=as.vector(sapply(time_points, function(x) rep(x, length(HSPC_populations)))),
                                   HSPC=rep(HSPC_populations, length(time_points)),
                                   stringsAsFactors = FALSE, check.names = FALSE)
  row.names(possible_names_mat) <- possible_names
  
  ### keep the required files only
  f_list <- f_list[which(f_list %in% paste0(possible_names, "_regress.Robj"))]
  
  ### for each of the Robj file, run the pseudotime analysis and combine the R objects for later use
  Combined_Seurat_Obj2 <- NULL
  for(f in f_list) {
    
    ### load the Seurat object and save the object name
    tmp_env <- new.env()
    load(paste0(Robj2_path, f), tmp_env)
    obj_name <- ls(tmp_env)
    assign("Seurat_Obj", get(obj_name, envir = tmp_env))
    rm(tmp_env)
    gc()
    
    ### rownames in the meta.data should be in the same order as colnames in the counts
    Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
    
    ### active assay = "RNA"
    Seurat_Obj@active.assay <- "RNA"
    
    ### annotate the file name to the R object
    tissue_name <- strsplit(f, split = "_", fixed = TRUE)[[1]][1]
    Seurat_Obj@meta.data$Tissue <- tissue_name
    
    ### annotate the development
    Seurat_Obj@meta.data$Development <- possible_names_mat[tissue_name, "Time"]
    
    ### annotate the HSPC type
    Seurat_Obj@meta.data$HSPC <- possible_names_mat[tissue_name, "HSPC"]
    
    ### combine the Seurat R objects
    if(is.null(Combined_Seurat_Obj2)) {
      Combined_Seurat_Obj2 <- Seurat_Obj
    } else {
      Combined_Seurat_Obj2 <- merge(x = Combined_Seurat_Obj2, y = Seurat_Obj,
                                    project = Seurat_Obj@project.name)
    }
    
  }
  
  ### annotate HSPC types in the original Seurat object
  Combined_Seurat_Obj@meta.data$HSPC <- "Other_Heme"
  Combined_Seurat_Obj@meta.data$HSPC[which(Combined_Seurat_Obj@meta.data$Cell_Type == "Stroma")] <- "Stroma"
  HSPC_barcodes <- intersect(Combined_Seurat_Obj@meta.data$GexCellFull,
                             Combined_Seurat_Obj2@meta.data$GexCellFull)
  Combined_Seurat_Obj@meta.data[HSPC_barcodes, "HSPC"] <- Combined_Seurat_Obj2@meta.data[HSPC_barcodes, "HSPC"]
  
  ### run PCA
  Combined_Seurat_Obj <- FindVariableFeatures(Combined_Seurat_Obj)
  Combined_Seurat_Obj <- ScaleData(Combined_Seurat_Obj)
  Combined_Seurat_Obj <- RunPCA(Combined_Seurat_Obj, npcs = 15)
  Combined_Seurat_Obj <- RunUMAP(Combined_Seurat_Obj, dims = 1:15)
  
  ### draw a PCA
  DimPlot(Combined_Seurat_Obj, reduction = "pca", group.by = "HSPC", pt.size = 1.5) +
    labs(title = paste0("PCA_Combined_Tissue"))
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Combined_Seurat_Obj)), rownames(Combined_Seurat_Obj@meta.data)))
  
  ### save the combined file as RDS
  saveRDS(Combined_Seurat_Obj,
          file = paste0(outputDir, "Combined_Seurat_Obj.RDS"))
  
}

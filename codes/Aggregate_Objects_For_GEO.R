###
#   File name : Aggregate_Objects_For_GEO.R
#   Author    : Hyunjin Kim
#   Date      : Jun 24, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : The seurat objects that were used in this project are spreaded.
#               Combine them and the annotations as well to submit 'Processed Data' to GEO
#   
#   Instruction
#               1. Source("Aggregate_Objects_For_GEO.R")
#               2. Run the function "aggregate_objects" - specify the input path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Aggregate_Objects_For_GEO.R/Aggregate_Objects_For_GEO.R")
#               > aggregate_objects(inputDir="Z:/ResearchHome/ProjectSpace/thomagrp/JCC282_Hematopoiesis/common/20210622 Objects for GEO/",
#                                   outputDir="./data/")
###

aggregate_objects <- function(inputDir="Z:/ResearchHome/ProjectSpace/thomagrp/JCC282_Hematopoiesis/common/20210622 Objects for GEO/",
                              outputDir="./data/GEO/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(data.table, quietly = TRUE)) {
    install.packages("data.table")
    require(data.table, quietly = TRUE)
  }
  
  ### create the output directory
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  
  ### loads an RData file, and returns it
  loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  
  ### load the original & combined seurat object
  original_seurat_obj <- loadRData(paste0(inputDir, "JCC282_Aggreg2.Robj"))
  gc()
    
  ### list of the seurat objects
  f <- list.files(path = inputDir,
                  pattern = "Final.Robj$")
  seurat_obj_list <- vector("list", length(f))
  names(seurat_obj_list) <- f
  
  ### load the annotated objects
  for(file in f) {
    seurat_obj_list[[file]] <- loadRData(paste0(inputDir, file))
  }
  
  ### see if the number of total cells in the original seurat object and the sumed up number from the annotated ones are the same
  print(paste("Sumed up from the annotated ones:", sum(sapply(seurat_obj_list, function(x) nrow(x@meta.data)))))
  print(paste("Number from the original seurat one:", nrow(original_seurat_obj@meta.data)))
  
  ### check cell names
  annotated_cell_names <- unlist(sapply(seurat_obj_list, function(x) rownames(x@meta.data)))
  original_cell_names <- rownames(original_seurat_obj@meta.data)
  
  ### merge the annotated objects
  combined_seurat_object <- NULL
  for(file in f) {
    if(is.null(combined_seurat_object)) {
      combined_seurat_object <- seurat_obj_list[[file]]
    } else {
      combined_seurat_object <- merge(x = combined_seurat_object, y = seurat_obj_list[[file]])
    }
  }
  print(paste("Number from the combined seurat object:", nrow(combined_seurat_object@meta.data)))
  
  ### check whether the orders are the same
  print(identical(rownames(combined_seurat_object@meta.data), colnames(combined_seurat_object@assays$RNA@counts)))
  print(identical(names(Idents(object = combined_seurat_object)), rownames(combined_seurat_object@meta.data)))
  
  ### remove unnecessary columns in the meta.data
  meta.data <- combined_seurat_object@meta.data
  colnames(meta.data)
  meta.data <- meta.data[,c("nCount_RNA", "nFeature_RNA", "nCount_ADT", "nFeature_ADT",
                            "Library", "Type", "Time", "Info",
                            "percent.mito", "percent.AlphaBeta", "percent.GammaDelta", "percent.IG",
                            "S.Score", "G2M.Score", "Phase", "hscScore", "Final_Annotation")]
  
  ### print out the count matrix and the meta.data
  fwrite(data.frame(Gene=rownames(combined_seurat_object@assays$RNA@counts),
                    combined_seurat_object@assays$RNA@counts,
                    stringsAsFactors = FALSE, check.names = FALSE),
         file = paste0(outputDir, "scRNASeq_Raw_Count_Matrix.tsv"),
         sep = "\t", row.names = FALSE)
  fwrite(data.frame(Cell_Barcode=rownames(meta.data),
                    meta.data,
                    stringsAsFactors = FALSE, check.names = FALSE),
         file = paste0(outputDir, "scRNASeq_Samples_Metadata.tsv"),
         sep = "\t", row.names = FALSE)
  
}

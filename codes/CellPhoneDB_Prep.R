###
#   File name : CellPhoneDB_Prep.R
#   Author    : Hyunjin Kim
#   Date      : Apr 14, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Make input txt files for CellPhoneDB from Trent's Seurat objects
#   
#   Instruction
#               1. Source("CellPhoneDB_Prep.R")
#               2. Run the function "cpdb_prep" - specify the input directory
#               3. The results will be generated under the corresponding directory
#
#   Example
#               > source("The_directory_of_CellPhoneDB_Prep.R/CellPhoneDB_Prep.R")
#               > cpdb_prep(Robj_dir="Z:/ResearchHome/ProjectSpace/thomagrp/JCC282_Hematopoiesis/common/RNA-Magnet/20210412_Signaling_updated_p-values/")
###

cpdb_prep <- function(Robj_dir="Z:/ResearchHome/ProjectSpace/thomagrp/JCC282_Hematopoiesis/common/RNA-Magnet/20210412_Signaling_updated_p-values/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  
  ### loads an RData file, and returns it
  loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  
  ### get sub directories
  subdirs <- list.dirs(Robj_dir, recursive = FALSE)
  
  ### for each sub directory
  for(sd in subdirs) {
    
    ### get the second sub directories
    subdirs2 <- list.dirs(sd, recursive = FALSE)
    
    ### for each sub directory2
    for(sd2 in subdirs2) {
      ### R obj path
      f <- list.files(sd2, pattern = "Specificity.Robj$")[1]
      
      ### load the R object
      Seurat_Obj <- loadRData(paste0(sd2, "/", f))
      gc()
      
      ### rownames in the meta.data should be in the same order as colnames in the counts
      Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
      print(identical(rownames(Seurat_Obj@meta.data), colnames(Seurat_Obj@assays$RNA@counts)))
      
      ### active assay = "RNA"
      Seurat_Obj@active.assay <- "RNA"
      
      ### write out the input count data
      write.table(data.frame(Gene=rownames(Seurat_Obj@assays$RNA@counts),
                             Seurat_Obj@assays$RNA@counts,
                             stringsAsFactors = FALSE, check.names = FALSE),
                  file = paste0(sd2, "/", strsplit(f, split = "_", fixed = TRUE)[[1]][1], "_counts.txt"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      
      ### write out the input meta data
      write.table(data.frame(Cell=rownames(Seurat_Obj@meta.data),
                             Cell_Type=as.character(Seurat_Obj@meta.data$HSPC),
                             stringsAsFactors = FALSE, check.names = FALSE),
                  file = paste0(sd2, "/", strsplit(f, split = "_", fixed = TRUE)[[1]][1], "_meta.txt"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      
      ### alert if there are more than 2 classes in the meta
      if(length(unique(Seurat_Obj@meta.data$HSPC)) > 2) {
        writeLines(paste("ERROR:", f, "There are more than 2 classes in HSPC column."))
      }
    }
    
  }
  
}

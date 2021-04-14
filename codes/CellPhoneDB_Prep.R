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
  
  ### from mouse gene names to human gene names based on gene orthology
  convert_MG_to_HG <- function(mouse_gene_symbols, method=c("gene_symbol", "ensembl_id")) {
    if(!require(biomaRt, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("biomaRt")
      require(biomaRt, quietly = TRUE)
    }
    
    ### get mouse-human map data
    mouse_db = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    mouse_human_gene_map <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "hsapiens_homolog_associated_gene_name", "hsapiens_homolog_ensembl_gene"), mart = mouse_db)
    
    ### only keep the genes that
    mouse_human_gene_map <- mouse_human_gene_map[which(mouse_human_gene_map$external_gene_name %in% mouse_gene_symbols),]
    
    ### if there is duplicated mouse genes, use the first one
    dups_idx <- which(duplicated(mouse_human_gene_map$external_gene_name))
    if(length(dups_idx) > 0) {
      mouse_human_gene_map <- mouse_human_gene_map[-dups_idx,]
    }
    
    ### only keep the rows with human genes
    if(method[1] == "gene_symbol") {
      mouse_human_gene_map <- mouse_human_gene_map[which(mouse_human_gene_map$hsapiens_homolog_associated_gene_name != ""),]
    } else {
      mouse_human_gene_map <- mouse_human_gene_map[which(mouse_human_gene_map$hsapiens_homolog_ensembl_gene != ""),]
    }
    
    ### set row names
    rownames(mouse_human_gene_map) <- mouse_human_gene_map$external_gene_name
    
    ### convert the mouse gene symbols to human gene symbols
    ### if none, return ""
    if(method[1] == "gene_symbol") {
      result <- mouse_human_gene_map[mouse_gene_symbols,"hsapiens_homolog_associated_gene_name"]
    } else {
      result <- mouse_human_gene_map[mouse_gene_symbols,"hsapiens_homolog_ensembl_gene"]
    }
    names(result) <- mouse_gene_symbols
    
    return(result)
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
      
      ### convert mouse gene names to human gene names
      h_genes <- convert_MG_to_HG(rownames(Seurat_Obj@assays$RNA@counts), method = "ensembl_id")
      h_genes <- h_genes[!is.na(h_genes)]
      
      ### write out the input count data
      write.table(data.frame(Gene=h_genes,
                             data.frame(Seurat_Obj@assays$RNA@counts, stringsAsFactors = FALSE, check.names = FALSE)[names(h_genes),],
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
      
      ### also write out the pre-built CellphoneDB command
      writeLines(paste0("cellphonedb method statistical_analysis ",
                        "/project_space/thomagrp/JCC282_Hematopoiesis/common/RNA-Magnet/20210412_Signaling_updated_p-values/Adult/LTHSC/AdultStromavLTHSC_meta.txt /project_space/thomagrp/JCC282_Hematopoiesis/common/RNA-Magnet/20210412_Signaling_updated_p-values/Adult/LTHSC/AdultStromavLTHSC_counts.txt --output-path=/project_space/thomagrp/JCC282_Hematopoiesis/common/RNA-Magnet/20210412_Signaling_updated_p-values/Adult/LTHSC/"))
    }
    
  }
  
}

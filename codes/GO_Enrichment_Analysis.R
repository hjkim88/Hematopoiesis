###
#   File name : GO_Enrichment_Analysis.R
#   Author    : Hyunjin Kim
#   Date      : Aug 5, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : GO Analyses on the clusters within each R object.
#               I’d like to know how the transcriptional programs of these immunophenotypic populations
#               change across development.
#   
#   Instruction
#               1. Source("GO_Enrichment_Analysis.R")
#               2. Run the function "go_analysis_trent" - specify the input directory and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_GO_Enrichment_Analysis.R/GO_Enrichment_Analysis.R")
#               > go_analysis_trent(Robj_path="C:/Users/hkim8/SJ/HSPC Analyses/HSPC Subsets/",
#                                   outputDir="./results/")
###

go_analysis_trent <- function(Robj_path="C:/Users/hkim8/SJ/HSPC Analyses/HSPC Subsets/",
                              outputDir="./results/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  options(java.parameters = "-Xmx10240m")
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  
  ### get Robj file list
  f_list <- list.files(Robj_path, pattern = ".Robj$")
  
  ### for each of the Robj file, run the GO analysis
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
    
    ### Find markers for each cluster?
    ### Or FindAllMarkers()?
    
    
    
  }
  
  
  
  
}

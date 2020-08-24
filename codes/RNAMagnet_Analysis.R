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
#               > rna_magnet_trent(Seurat_RObj_path="./data/Combined_BM_Seurat_Obj.RDATA",
#                                  outputDir="./results/RNA-Magnet/")
###

rna_magnet_trent <- function(Seurat_RObj_path="./data/Combined_BM_Seurat_Obj.RDATA",
                             outputDir="./results/RNA-Magnet/") {
  
  ### load libraries
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
  
  ### load the Seurat object and save the object name
  tmp_env <- new.env()
  load(Seurat_RObj_path, tmp_env)
  obj_name <- ls(tmp_env)
  assign("Seurat_Obj", get(obj_name, envir = tmp_env))
  rm(tmp_env)
  gc()
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data)))
  
  
  
  
  


}

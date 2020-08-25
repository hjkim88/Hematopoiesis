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
  if(!require(reticulate, quietly = TRUE)) {
    install.packages("reticulate")
    require(reticulate, quietly = TRUE)
  }
  if(!require(Rmagic, quietly = TRUE)) {
    install.packages("Rmagic")
    require(Rmagic, quietly = TRUE)
  }
  
  ### see python environment since RNAMagnet uses the python module 'magic'
  # conda_create("r-reticulate")
  # conda_install("r-reticulate", "python-magic")
  use_condaenv("r-reticulate")
  py_config()
  py_module_available("magic")
  
  ### load the Seurat object and save the object name
  tmp_env <- new.env()
  load(Seurat_RObj_path, tmp_env)
  obj_name <- ls(tmp_env)
  assign("Seurat_Obj", get(obj_name, envir = tmp_env))
  rm(tmp_env)
  gc()
  
  ### determine necessary variables
  time_points <- c("E16", "E18", "P0", "ADULT")
  cell_type <- c("Heme", "Stroma")
  
  ### active assay = "RNA"
  Seurat_Obj@active.assay <- "RNA"
  
  ### preprocessing for PCA and UMAP
  Seurat_Obj <- FindVariableFeatures(Seurat_Obj)
  Seurat_Obj <- ScaleData(Seurat_Obj)
  
  ### run PCA
  Seurat_Obj <- RunPCA(Seurat_Obj, npcs = 15)
  
  ### run UMAP
  Seurat_Obj <- RunUMAP(Seurat_Obj, dims = 1:5)
  
  ### order the data by time
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[order(factor(Seurat_Obj@meta.data$Development,
                                                            levels = time_points)),]
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
  
  ### set the ident of the object with the HSPC type
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$Tissue)
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data)))
  
  ### for each R object, perform RNAMagnet
  for(obj in unique(Seurat_Obj@meta.data$Tissue)) {
    
    ### new output directory
    outputDir2 <- paste0(outputDir, obj, "/")
    dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
    
    ### split the Seurat obj based on obj info
    subset_Seurat_Obj <- subset(Seurat_Obj, idents=obj)
    
    ### check whether the orders are the same
    print(identical(names(Idents(object = subset_Seurat_Obj)), rownames(subset_Seurat_Obj@meta.data)))
    
    ### rownames in the meta.data should be in the same order as colnames in the counts
    subset_Seurat_Obj@meta.data <- subset_Seurat_Obj@meta.data[colnames(subset_Seurat_Obj@assays$RNA@counts),]
    
    ### There should be at least two clusters for RNAMagnet
    if(length(unique(subset_Seurat_Obj@meta.data$seurat_clusters)) > 1) {
      ### for each cluster, set the cluster as the anchor and run RNAMagnet
      for(clust in unique(subset_Seurat_Obj@meta.data$seurat_clusters)) {
        
        ### set the ident of the seurat object with the cluster info
        subset_Seurat_Obj <- SetIdent(object = subset_Seurat_Obj,
                                      cells = rownames(subset_Seurat_Obj@meta.data),
                                      value = subset_Seurat_Obj@meta.data$seurat_clusters)
        
        ### run RNAMagnet
        result <- RNAMagnetAnchors(subset_Seurat_Obj,
                                   anchors = unique(subset_Seurat_Obj@meta.data$seurat_clusters))
        
        ### RNAMagnet on PCA
        qplot(x =Embeddings(subset_Seurat_Obj,reduction="pca")[,1],
              y=Embeddings(subset_Seurat_Obj,reduction="pca")[,2],
              color = direction,
              size=I(0.75),
              alpha= adhesiveness,data=result) +
          scale_color_brewer(name = "RNAMagnet\nLocation",palette= "Set1") +
          scale_alpha_continuous(name = "RNAMagnet\nAdhesiveness") +
          theme_bw() +
          theme(panel.grid = element_blank(),
                axis.text = element_blank(),
                axis.title = element_blank())
        
      }
    }
    
  }
  
}

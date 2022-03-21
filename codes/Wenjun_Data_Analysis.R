###
#   File name : Wenjun_Data_Analysis.R
#   Author    : Hyunjin Kim
#   Date      : Mar 18, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Cell annotation first and run RNAMagnet/CellComm to look at the interactions among different cell types
#
#   CAR (CXCL12-Abundant Reticular) cell marker: CXCL12, SCF, FOXC1, EBF3
#
#   Instruction
#               1. Source("Wenjun_Data_Analysis.R")
#               2. Run the function "wenjun_preprocess" - specify the input path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Wenjun_Data_Analysis.R/Wenjun_Data_Analysis.R")
#               > wenjun_preprocess(Seurat_Obj_Path="./data/Wenjun_Seurat_Obj.RDS",
#                                   outputDir="./results/Wenjun_Results/")
###

wenjun_analysis <- function(Seurat_Obj_Path="./data/Wenjun_Seurat_Obj.RDS",
                            outputDir="./results/Wenjun_Results/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(SingleR, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("SingleR")
    require(SingleR, quietly = TRUE)
  }
  if(!require(celldex, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("celldex")
    require(celldex, quietly = TRUE)
  }
  if(!require(scran, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("scran")
    require(scran, quietly = TRUE)
  }
  if(!require(RNAMagnet, quietly = TRUE)) {
    remotes::install_github("veltenlab/rnamagnet")
    require(RNAMagnet, quietly = TRUE)
  }
  if(!require(patchwork, quietly = TRUE)) {
    install.packages("patchwork")
    require(patchwork, quietly = TRUE)
  }
  if(!require(NMF, quietly = TRUE)) {
    install.packages("NMF")
    require(NMF, quietly = TRUE)
  }
  if(!require(circlize, quietly = TRUE)) {
    devtools::install_github("jokergoo/circlize")
    require(circlize, quietly = TRUE)
  }
  if(!require(ComplexHeatmap, quietly = TRUE)) {
    devtools::install_github("jokergoo/ComplexHeatmap")
    require(ComplexHeatmap, quietly = TRUE)
  }
  if(!require(CellChat, quietly = TRUE)) {
    devtools::install_github("sqjin/CellChat")
    require(CellChat, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  
  ### load seurat object
  WJ_Seurat_Obj <- readRDS(file = Seurat_Obj_Path)
  
  ### check whether the orders are the same
  print(identical(rownames(WJ_Seurat_Obj@meta.data), colnames(WJ_Seurat_Obj@assays$RNA@counts)))
  print(identical(names(Idents(object = WJ_Seurat_Obj)), rownames(WJ_Seurat_Obj@meta.data)))
  
  ### set markers to use for the annotation
  # ref <- MouseRNAseqData()
  # ref <- ImmGenData()
  ref <- BlueprintEncodeData()
  markers <- scran::findMarkers(ref, groups = ref$label.main, 
                                test.type = "t", assay.type = "logcounts",
                                lfc = 0.6, pval.type = "all", direction = "up")
  markerlist <- lapply(markers, function(w) rownames(w)[w$FDR < 0.05])
  sapply(markerlist, length)
  
  ### annotation using singleR
  pred.sce <- SingleR::SingleR(test = WJ_Seurat_Obj@assays$RNA@data, ref = ref, labels = ref$label.main)
  WJ_Seurat_Obj$singleR_celltypes_Blueprint <- pred.sce$pruned.labels
  cell_num <- sapply(unique(WJ_Seurat_Obj$singleR_celltypes_Blueprint), function(x) {
    return(length(which(WJ_Seurat_Obj$singleR_celltypes_Blueprint == x)))
  })
  
  ### cells that have at least 300 cells
  target_types <- unlist(sapply(names(cell_num), function(x) {
    if(!is.na(x) && cell_num[x] > 300) {
      return(x)
    } else {
      return(NULL)
    }
  }, USE.NAMES = FALSE))
  WJ_Seurat_Obj$singleR_celltypes_Blueprint2 <- WJ_Seurat_Obj$singleR_celltypes_Blueprint
  WJ_Seurat_Obj$singleR_celltypes_Blueprint2[which(!WJ_Seurat_Obj$singleR_celltypes_Blueprint %in% target_types)] <- NA
  
  
  DimPlot(object = WJ_Seurat_Obj, reduction = "umap", raster = TRUE,
          group.by = "Cell1",
          pt.size = 1)
  DimPlot(object = WJ_Seurat_Obj, reduction = "umap", raster = TRUE,
          group.by = "Cell2",
          pt.size = 1)
  
  
  ### UMAP of the annotations
  p <- DimPlot(object = WJ_Seurat_Obj, reduction = "umap", raster = TRUE,
          group.by = "singleR_celltypes_Blueprint",
          pt.size = 1) +
    ggtitle(paste0("Annotated Cell Types")) +
    labs(color = "Annotation") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 48) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 36, color = "black", face = "bold"),
          axis.title = element_text(size = 36, color = "black", face = "bold"),
          axis.text = element_text(size = 36, color = "black", face = "bold"),
          legend.title = element_text(size = 24, color = "black", face = "bold"),
          legend.text = element_text(size = 24, color = "black", face = "bold"),
          axis.ticks = element_blank())
  p <- DimPlot(object = WJ_Seurat_Obj, reduction = "umap", raster = TRUE,
               group.by = "singleR_celltypes_Blueprint2",
               pt.size = 1) +
    ggtitle(paste0("Annotated Cell Types")) +
    labs(color = "Annotation") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 48) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 36, color = "black", face = "bold"),
          axis.title = element_text(size = 36, color = "black", face = "bold"),
          axis.text = element_text(size = 36, color = "black", face = "bold"),
          legend.title = element_text(size = 24, color = "black", face = "bold"),
          legend.text = element_text(size = 24, color = "black", face = "bold"),
          axis.ticks = element_blank())
  
  ### CXCL12, FOXC1, EBF3 UMAP
  p <- FeaturePlot(WJ_Seurat_Obj, features = c("Cxcl12", "Foxc1", "Ebf3"), cols = c("lightgray", "red"))
  
  ### cluster approach
  DimPlot(object = WJ_Seurat_Obj, reduction = "umap", raster = TRUE,
          group.by = "seurat_clusters",
          pt.size = 1)
  WJ_Seurat_Obj <- SetIdent(object = WJ_Seurat_Obj,
                            cells = rownames(WJ_Seurat_Obj@meta.data),
                            value = WJ_Seurat_Obj$seurat_clusters)
  de_result_all <- FindAllMarkers(WJ_Seurat_Obj,
                                  min.pct = 0.2,
                                  logfc.threshold = 0.2,
                                  test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result_all),
                         de_result_all,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/Wenjun_Cluster_AllMarkers.xlsx"),
              sheetName = "Wenjun_Cluster_AllMarkers", row.names = FALSE)
  
  p <- DimPlot(object = WJ_Seurat_Obj, reduction = "umap", raster = TRUE,
               group.by = "seurat_clusters",
               pt.size = 0.5, label = TRUE, label.size = 12) +
    ggtitle(paste0("")) +
    labs(color = "Clusters") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 48) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 36, color = "black", face = "bold"),
          axis.title = element_text(size = 36, color = "black", face = "bold"),
          axis.text = element_text(size = 36, color = "black", face = "bold"),
          legend.title = element_text(size = 24, color = "black", face = "bold"),
          legend.text = element_text(size = 24, color = "black", face = "bold"),
          axis.ticks = element_blank())
  ggsave(paste0(outputDir, "/Wenjun_UMAP_Clusters_label.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  
  ### Cellchat
  
  ### create a cell chat object
  cellchat <- createCellChat(object = WJ_Seurat_Obj, group.by = "groups")
  
  ### set ligand-receptor interaction db
  # CellChatDB <- CellChatDB.human
  CellChatDB <- CellChatDB.mouse
  
  # use a subset of CellChatDB for cell-cell communication analysis
  # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  # use all CellChatDB for cell-cell communication analysis
  CellChatDB.use <- CellChatDB
  
  ### set the used database in the object
  cellchat@DB <- CellChatDB.use
  
  ### subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  
  ### Preprocessing
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # project gene expression data onto PPI network (optional)
  cellchat <- projectData(cellchat, PPI.human)
  
  ### Compute the communication probability and infer cellular communication network
  cellchat <- computeCommunProb(cellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  ### Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  
  ### Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  
  
}

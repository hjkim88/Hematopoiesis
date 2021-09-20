###
#   File name : Additional_Analysis_Sep_2021.R
#   Author    : Hyunjin Kim
#   Date      : Sep 20, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Look at the TDH027 - p21 stroma data - sub-cluster
#               a. If there is a cluster that has some genes and UMIs that are comparable to other seurat object (e.g., SJCAR19), then we are fine
#               b. if they appear to be cycling less or if a cluster has less genes/cell is cycling less than other clusters
#                  (cycling less would be in G1 or G0) so investigate the correlation between cycling (cycling gene expression or phase)
#                  and the number of genes/cells
#   
#   Instruction
#               1. Source("Additional_Analysis_Sep_2021.R")
#               2. Run the function "sub_cluster_checking" - specify the input path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Additional_Analysis_Sep_2021.R/Additional_Analysis_Sep_2021.R")
#               > sub_cluster_checking(CellRangerPath="Z:/ResearchHome/SharedResources/Immunoinformatics/common/JCC/JCC395/JCC392_TDH027/filtered_feature_bc_matrix/")
###

sub_cluster_checking <- function(CellRangerPath="Z:/ResearchHome/SharedResources/Immunoinformatics/common/JCC/JCC395/JCC392_TDH027/filtered_feature_bc_matrix/",
                                 ExampleDataPath="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(dplyr, quietly = TRUE)) {
    install.packages("dplyr")
    require(dplyr, quietly = TRUE)
  }
  
  ### load the TDH027 dataset
  TDH027.data <- Read10X(data.dir = CellRangerPath)
  
  ### create a Seurat object
  ### min.cells : include genes detected in at least this many cells
  ### min.features : include cells where at least this many features are detected
  TDH027_Seurat_Obj <- CreateSeuratObject(counts = TDH027.data$`Gene Expression`,
                                          project = "SMF_TDH027",
                                          min.cells = 3,
                                          min.features = 200)
  
  ### add ADT assay
  TDH027_Seurat_Obj[["ADT"]] <- CreateAssayObject(counts = TDH027.data$`Antibody Capture`[, colnames(TDH027_Seurat_Obj)])
  
  ### active assay = "RNA"
  TDH027_Seurat_Obj@active.assay <- "RNA"
  
  ### MT percentage
  TDH027_Seurat_Obj[["percent.mt"]] <- PercentageFeatureSet(TDH027_Seurat_Obj, pattern = "^MT-")
  
  ### Visualize QC metrics as a violin plot
  VlnPlot(TDH027_Seurat_Obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot(density(TDH027_Seurat_Obj@meta.data$nFeature_RNA))
  plot(density(TDH027_Seurat_Obj@meta.data$percent.mt))
  
  ### Filter cells that have unique feature counts over 6500 or less than 300
  ### Filter cells that have > 10% mitochondrial counts
  TDH027_Seurat_Obj <- subset(TDH027_Seurat_Obj, subset = nFeature_RNA > 300 & nFeature_RNA < 6500 & percent.mt < 10)
  
  ### Cell cycle score (will be used later for regression out)
  TDH027_Seurat_Obj <- CellCycleScoring(object = TDH027_Seurat_Obj,
                                        g2m.features = cc.genes$g2m.genes,
                                        s.features = cc.genes$s.genes)
  ### normalization
  TDH027_Seurat_Obj <- NormalizeData(TDH027_Seurat_Obj,
                                     normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  TDH027_Seurat_Obj <- FindVariableFeatures(TDH027_Seurat_Obj,
                                            selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  TDH027_Seurat_Obj <- ScaleData(TDH027_Seurat_Obj,
                                 vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### PCA
  TDH027_Seurat_Obj <- RunPCA(TDH027_Seurat_Obj,
                              features = VariableFeatures(object = TDH027_Seurat_Obj))
  
  ### UMAP
  TDH027_Seurat_Obj <- RunUMAP(TDH027_Seurat_Obj, dims = 1:15)
  
  ### clustering
  TDH027_Seurat_Obj <- FindNeighbors(TDH027_Seurat_Obj, dims = 1:10)
  TDH027_Seurat_Obj <- FindClusters(TDH027_Seurat_Obj, resolution = 0.5)
  
  ### check whether the orders are the same
  print(identical(rownames(TDH027_Seurat_Obj@meta.data), colnames(TDH027_Seurat_Obj@assays$RNA@counts)))
  
  ### check nFeature_RNA
  print(TDH027_Seurat_Obj$nFeature_RNA[1])
  print(length(which(TDH027_Seurat_Obj@assays$RNA@counts[,1] > 0)))
  
  ### see the UMAP plot
  DimPlot(object = TDH027_Seurat_Obj, reduction = "umap",
          group.by = "seurat_clusters", label = TRUE,
          pt.size = 0.5)
  FeaturePlot(TDH027_Seurat_Obj, features = c("Cxcl12"), cols = c("lightgray", "red"), label = TRUE, pt.size = 1.5)
  FeaturePlot(TDH027_Seurat_Obj, features = c("nFeature_RNA"), cols = c("lightgray", "red"), label = TRUE, pt.size = 1.5)
  
  ### which clusters have less genes/UMIs?
  clstr_gene_num <- rep(0, length(levels(TDH027_Seurat_Obj$seurat_clusters)))
  names(clstr_gene_num) <- levels(TDH027_Seurat_Obj$seurat_clusters)
  
  for(clstr in levels(TDH027_Seurat_Obj$seurat_clusters)) {
    clstr_gene_num[clstr] <- mean(TDH027_Seurat_Obj$nFeature_RNA[which(TDH027_Seurat_Obj$seurat_clusters == clstr)])
  }
  
  ### load Jeremy's object
  JCC_Seurat_Obj <- readRDS(file = ExampleDataPath)
  
  
  
  
}

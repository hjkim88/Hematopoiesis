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
#               > sub_cluster_checking(CellRangerPath1="Z:/ResearchHome/SharedResources/Immunoinformatics/common/JCC/JCC395/JCC392_TDH027/filtered_feature_bc_matrix/")
###

sub_cluster_checking <- function(CellRangerPath1="Z:/ResearchHome/SharedResources/Immunoinformatics/common/JCC/JCC395/JCC392_TDH027/filtered_feature_bc_matrix/",
                                 CellRangerPath2="Z:/ResearchHome/SharedResources/Immunoinformatics/common/JCC/JCC395/JCC392_TDH019/filtered_feature_bc_matrix/",
                                 ExampleDataPath="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds",
                                 outputDir="./results/Sep_2021/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(dplyr, quietly = TRUE)) {
    install.packages("dplyr")
    require(dplyr, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(ggrepel, quietly = TRUE)) {
    install.packages("ggrepel")
    require(ggrepel, quietly = TRUE)
  }
  
  ### load the TDH027 dataset
  TDH027.data <- Read10X(data.dir = CellRangerPath1)
  
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
  
  
  ### load the TDH019 dataset
  TDH019.data <- Read10X(data.dir = CellRangerPath2)
  
  ### create a Seurat object
  ### min.cells : include genes detected in at least this many cells
  ### min.features : include cells where at least this many features are detected
  TDH019_Seurat_Obj <- CreateSeuratObject(counts = TDH019.data$`Gene Expression`,
                                          project = "SMF_TDH019",
                                          min.cells = 3,
                                          min.features = 200)
  
  ### add ADT assay
  TDH019_Seurat_Obj[["ADT"]] <- CreateAssayObject(counts = TDH019.data$`Antibody Capture`[, colnames(TDH019_Seurat_Obj)])
  
  ### active assay = "RNA"
  TDH019_Seurat_Obj@active.assay <- "RNA"
  
  ### MT percentage
  TDH019_Seurat_Obj[["percent.mt"]] <- PercentageFeatureSet(TDH019_Seurat_Obj, pattern = "^MT-")
  
  ### Visualize QC metrics as a violin plot
  VlnPlot(TDH019_Seurat_Obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot(density(TDH019_Seurat_Obj@meta.data$nFeature_RNA))
  plot(density(TDH019_Seurat_Obj@meta.data$percent.mt))
  
  ### Filter cells that have unique feature counts over 6500 or less than 300
  ### Filter cells that have > 10% mitochondrial counts
  TDH019_Seurat_Obj <- subset(TDH019_Seurat_Obj, subset = nFeature_RNA > 300 & nFeature_RNA < 6500 & percent.mt < 10)
  
  ### Cell cycle score (will be used later for regression out)
  TDH019_Seurat_Obj <- CellCycleScoring(object = TDH019_Seurat_Obj,
                                        g2m.features = cc.genes$g2m.genes,
                                        s.features = cc.genes$s.genes)
  ### normalization
  TDH019_Seurat_Obj <- NormalizeData(TDH019_Seurat_Obj,
                                     normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  TDH019_Seurat_Obj <- FindVariableFeatures(TDH019_Seurat_Obj,
                                            selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  TDH019_Seurat_Obj <- ScaleData(TDH019_Seurat_Obj,
                                 vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### PCA
  TDH019_Seurat_Obj <- RunPCA(TDH019_Seurat_Obj,
                              features = VariableFeatures(object = TDH019_Seurat_Obj))
  
  ### UMAP
  TDH019_Seurat_Obj <- RunUMAP(TDH019_Seurat_Obj, dims = 1:15)
  
  ### clustering
  TDH019_Seurat_Obj <- FindNeighbors(TDH019_Seurat_Obj, dims = 1:10)
  TDH019_Seurat_Obj <- FindClusters(TDH019_Seurat_Obj, resolution = 0.5)
  
  ### check whether the orders are the same
  print(identical(rownames(TDH019_Seurat_Obj@meta.data), colnames(TDH019_Seurat_Obj@assays$RNA@counts)))
  
  ### check nFeature_RNA
  print(TDH019_Seurat_Obj$nFeature_RNA[1])
  print(length(which(TDH019_Seurat_Obj@assays$RNA@counts[,1] > 0)))
  
  
  ### cell cycle color
  cell_cycle_color <- c("lightgray", "orange", "red")
  names(cell_cycle_color) <- c("G1", "S", "G2M")
  
  
  ### create the output directory
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  
  ### see the UMAP plot
  p <- DimPlot(object = TDH027_Seurat_Obj, reduction = "umap",
          group.by = "seurat_clusters", label = TRUE,
          pt.size = 1.5, label.size = 10) +
    ggtitle("TDH027 Clusters") +
    theme_classic(base_size = 50) +
    theme(legend.position = "none")
  ggsave(paste0(outputDir, "TDH027_Clusters_UMAP.png"), plot = p, width = 15, height = 12, dpi = 350)
  
  p <- FeaturePlot(TDH027_Seurat_Obj, features = c("Cxcl12"),
                   cols = c("lightgray", "red"), label = TRUE, pt.size = 1.5, label.size = 10) +
    ggtitle("TDH027 Cxcl12 Expression") +
    theme_classic(base_size = 50)
  ggsave(paste0(outputDir, "TDH027_Cxcl12_UMAP.png"), plot = p, width = 15, height = 12, dpi = 350)
  
  p <- FeaturePlot(TDH027_Seurat_Obj, features = c("Cxcr4"),
                   cols = c("lightgray", "red"), label = TRUE, pt.size = 1.5, label.size = 10) +
    ggtitle("TDH027 Cxcr4 Expression") +
    theme_classic(base_size = 50)
  ggsave(paste0(outputDir, "TDH027_Cxcr4_UMAP.png"), plot = p, width = 15, height = 12, dpi = 350)
  
  p <- FeaturePlot(TDH027_Seurat_Obj, features = c("nFeature_RNA"),
                   cols = c("lightgray", "red"), label = TRUE, pt.size = 1.5, label.size = 10) +
    ggtitle("TDH027 Gene #") +
    theme_classic(base_size = 50)
  ggsave(paste0(outputDir, "TDH027_nFeature_UMAP.png"), plot = p, width = 15, height = 12, dpi = 350)
  
  p <- FeaturePlot(TDH027_Seurat_Obj, features = c("Phase"),
                   label = TRUE, pt.size = 1.5, label.size = 10) +
    ggtitle("TDH027 Cell Cycle") +
    scale_color_manual(values = cell_cycle_color) +
    theme_classic(base_size = 50)
  ggsave(paste0(outputDir, "TDH027_Cell_Cycle_UMAP.png"), plot = p, width = 15, height = 12, dpi = 350)
  
  ### which clusters have less genes/UMIs?
  TDH027_clstr_gene_num <- rep(0, length(levels(TDH027_Seurat_Obj$seurat_clusters)))
  names(TDH027_clstr_gene_num) <- levels(TDH027_Seurat_Obj$seurat_clusters)
  TDH027_clstr_s_score <- rep(0, length(levels(TDH027_Seurat_Obj$seurat_clusters)))
  names(TDH027_clstr_s_score) <- levels(TDH027_Seurat_Obj$seurat_clusters)
  TDH027_clstr_g2m_score <- rep(0, length(levels(TDH027_Seurat_Obj$seurat_clusters)))
  names(TDH027_clstr_g2m_score) <- levels(TDH027_Seurat_Obj$seurat_clusters)
  
  for(clstr in levels(TDH027_Seurat_Obj$seurat_clusters)) {
    TDH027_clstr_gene_num[clstr] <- mean(TDH027_Seurat_Obj$nFeature_RNA[which(TDH027_Seurat_Obj$seurat_clusters == clstr)])
    TDH027_clstr_s_score[clstr] <- mean(TDH027_Seurat_Obj$S.Score[which(TDH027_Seurat_Obj$seurat_clusters == clstr)])
    TDH027_clstr_g2m_score[clstr] <- mean(TDH027_Seurat_Obj$G2M.Score[which(TDH027_Seurat_Obj$seurat_clusters == clstr)])
  }
  
  
  ### TDH019
  p <- DimPlot(object = TDH019_Seurat_Obj, reduction = "umap",
               group.by = "seurat_clusters", label = TRUE,
               pt.size = 1.5, label.size = 10) +
    ggtitle("TDH019 Clusters") +
    theme_classic(base_size = 50) +
    theme(legend.position = "none")
  ggsave(paste0(outputDir, "TDH019_Clusters_UMAP.png"), plot = p, width = 15, height = 12, dpi = 350)
  
  p <- FeaturePlot(TDH019_Seurat_Obj, features = c("Cxcl12"),
                   cols = c("lightgray", "red"), label = TRUE, pt.size = 1.5, label.size = 10) +
    ggtitle("TDH019 Cxcl12 Expression") +
    theme_classic(base_size = 50)
  ggsave(paste0(outputDir, "TDH019_Cxcl12_UMAP.png"), plot = p, width = 15, height = 12, dpi = 350)
  
  p <- FeaturePlot(TDH019_Seurat_Obj, features = c("Cxcr4"),
                   cols = c("lightgray", "red"), label = TRUE, pt.size = 1.5, label.size = 10) +
    ggtitle("TDH019 Cxcr4 Expression") +
    theme_classic(base_size = 50)
  ggsave(paste0(outputDir, "TDH019_Cxcr4_UMAP.png"), plot = p, width = 15, height = 12, dpi = 350)
  
  p <- FeaturePlot(TDH019_Seurat_Obj, features = c("nFeature_RNA"),
                   cols = c("lightgray", "red"), label = TRUE, pt.size = 1.5, label.size = 10) +
    ggtitle("TDH019 Gene #") +
    theme_classic(base_size = 50)
  ggsave(paste0(outputDir, "TDH019_nFeature_UMAP.png"), plot = p, width = 15, height = 12, dpi = 350)
  
  p <- FeaturePlot(TDH019_Seurat_Obj, features = c("Phase"),
                   label = TRUE, pt.size = 1.5, label.size = 10) +
    ggtitle("TDH019 Cell Cycle") +
    scale_color_manual(values = cell_cycle_color) +
    theme_classic(base_size = 50)
  ggsave(paste0(outputDir, "TDH019_Cell_Cycle_UMAP.png"), plot = p, width = 15, height = 12, dpi = 350)
  
  ### which clusters have less genes/UMIs?
  TDH019_clstr_gene_num <- rep(0, length(levels(TDH019_Seurat_Obj$seurat_clusters)))
  names(TDH019_clstr_gene_num) <- levels(TDH019_Seurat_Obj$seurat_clusters)
  TDH019_clstr_s_score <- rep(0, length(levels(TDH019_Seurat_Obj$seurat_clusters)))
  names(TDH019_clstr_s_score) <- levels(TDH019_Seurat_Obj$seurat_clusters)
  TDH019_clstr_g2m_score <- rep(0, length(levels(TDH019_Seurat_Obj$seurat_clusters)))
  names(TDH019_clstr_g2m_score) <- levels(TDH019_Seurat_Obj$seurat_clusters)
  
  for(clstr in levels(TDH019_Seurat_Obj$seurat_clusters)) {
    TDH019_clstr_gene_num[clstr] <- mean(TDH019_Seurat_Obj$nFeature_RNA[which(TDH019_Seurat_Obj$seurat_clusters == clstr)])
    TDH019_clstr_s_score[clstr] <- mean(TDH019_Seurat_Obj$S.Score[which(TDH019_Seurat_Obj$seurat_clusters == clstr)])
    TDH019_clstr_g2m_score[clstr] <- mean(TDH019_Seurat_Obj$G2M.Score[which(TDH019_Seurat_Obj$seurat_clusters == clstr)])
  }
  
  
  ### correlation between nFeature_RNA vs Cell cycle score
  plot_df <- data.frame(Cluster=names(TDH027_clstr_gene_num),
                        TDH027_Gene_Num=TDH027_clstr_gene_num,
                        TDH027_S_Score=TDH027_clstr_s_score,
                        TDH027_G2M_Score=TDH027_clstr_g2m_score,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  p_cor <- round(cor(plot_df$TDH027_Gene_Num,
                     plot_df$TDH027_S_Score, method = "spearman", use = "complete.obs"), 2)
  pv <- round(cor.test(plot_df$TDH027_Gene_Num,
                       plot_df$TDH027_S_Score, method = "spearman", use = "complete.obs")$p.value, 2)
  p <- ggplot(data = plot_df, aes(x=TDH027_Gene_Num, y=TDH027_S_Score)) +
    geom_point(col = "#487A8F", size = 8) +
    labs(title = paste0("Spearman Correlation:", p_cor),
         subtitle = paste0("P-value:", pv)) +
    xlab("TDH027 Average Gene #") +
    ylab("TDH027 S Score") +
    geom_label_repel(aes(label = Cluster),
                     size = 5,
                     col = "#3B3B53",
                     segment.color = "#3B3B53",
                     nudge_x = 0.03) +
    geom_smooth(method = lm, color="#AA4C26", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 40))
  ggsave(filename = paste0(outputDir, "TDH027_Correlation_Gene_Num_S_Score.png"), plot = p, width = 12, height = 10, dpi = 350)
  
  p_cor <- round(cor(plot_df$TDH027_Gene_Num,
                     plot_df$TDH027_G2M_Score, method = "spearman", use = "complete.obs"), 2)
  pv <- round(cor.test(plot_df$TDH027_Gene_Num,
                       plot_df$TDH027_G2M_Score, method = "spearman", use = "complete.obs")$p.value, 2)
  p <- ggplot(data = plot_df, aes(x=TDH027_Gene_Num, y=TDH027_G2M_Score)) +
    geom_point(col = "#487A8F", size = 8) +
    labs(title = paste0("Spearman Correlation:", p_cor),
         subtitle = paste0("P-value:", pv)) +
    xlab("TDH027 Average Gene #") +
    ylab("TDH027 G2M Score") +
    geom_label_repel(aes(label = Cluster),
                     size = 5,
                     col = "#3B3B53",
                     segment.color = "#3B3B53",
                     nudge_x = 0.03) +
    geom_smooth(method = lm, color="#AA4C26", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 40))
  ggsave(filename = paste0(outputDir, "TDH027_Correlation_Gene_Num_G2M_Score.png"), plot = p, width = 12, height = 10, dpi = 350)
  
  
  ### TDH019
  plot_df <- data.frame(Cluster=names(TDH019_clstr_gene_num),
                        TDH019_Gene_Num=TDH019_clstr_gene_num,
                        TDH019_S_Score=TDH019_clstr_s_score,
                        TDH019_G2M_Score=TDH019_clstr_g2m_score,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  p_cor <- round(cor(plot_df$TDH019_Gene_Num,
                     plot_df$TDH019_S_Score, method = "spearman", use = "complete.obs"), 2)
  pv <- round(cor.test(plot_df$TDH019_Gene_Num,
                       plot_df$TDH019_S_Score, method = "spearman", use = "complete.obs")$p.value, 2)
  p <- ggplot(data = plot_df, aes(x=TDH019_Gene_Num, y=TDH019_S_Score)) +
    geom_point(col = "#487A8F", size = 8) +
    labs(title = paste0("Spearman Correlation:", p_cor),
         subtitle = paste0("P-value:", pv)) +
    xlab("TDH019 Average Gene #") +
    ylab("TDH019 S Score") +
    geom_label_repel(aes(label = Cluster),
                     size = 5,
                     col = "#3B3B53",
                     segment.color = "#3B3B53",
                     nudge_x = 0.03) +
    geom_smooth(method = lm, color="#AA4C26", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 40))
  ggsave(filename = paste0(outputDir, "TDH019_Correlation_Gene_Num_S_Score.png"), plot = p, width = 12, height = 10, dpi = 350)
  
  p_cor <- round(cor(plot_df$TDH019_Gene_Num,
                     plot_df$TDH019_G2M_Score, method = "spearman", use = "complete.obs"), 2)
  pv <- round(cor.test(plot_df$TDH019_Gene_Num,
                       plot_df$TDH019_G2M_Score, method = "spearman", use = "complete.obs")$p.value, 2)
  p <- ggplot(data = plot_df, aes(x=TDH019_Gene_Num, y=TDH019_G2M_Score)) +
    geom_point(col = "#487A8F", size = 8) +
    labs(title = paste0("Spearman Correlation:", p_cor),
         subtitle = paste0("P-value:", pv)) +
    xlab("TDH019 Average Gene #") +
    ylab("TDH019 G2M Score") +
    geom_label_repel(aes(label = Cluster),
                     size = 5,
                     col = "#3B3B53",
                     segment.color = "#3B3B53",
                     nudge_x = 0.03) +
    geom_smooth(method = lm, color="#AA4C26", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 40))
  ggsave(filename = paste0(outputDir, "TDH019_Correlation_Gene_Num_G2M_Score.png"), plot = p, width = 12, height = 10, dpi = 350)
  
  #
  ### removing the smallest cluster
  #
  
  ### correlation between nFeature_RNA vs Cell cycle score
  plot_df <- data.frame(Cluster=names(TDH027_clstr_gene_num)[-13],
                        TDH027_Gene_Num=TDH027_clstr_gene_num[-13],
                        TDH027_S_Score=TDH027_clstr_s_score[-13],
                        TDH027_G2M_Score=TDH027_clstr_g2m_score[-13],
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  p_cor <- round(cor(plot_df$TDH027_Gene_Num,
                     plot_df$TDH027_S_Score, method = "spearman", use = "complete.obs"), 2)
  pv <- round(cor.test(plot_df$TDH027_Gene_Num,
                       plot_df$TDH027_S_Score, method = "spearman", use = "complete.obs")$p.value, 2)
  p <- ggplot(data = plot_df, aes(x=TDH027_Gene_Num, y=TDH027_S_Score)) +
    geom_point(col = "#487A8F", size = 8) +
    labs(title = paste0("Spearman Correlation:", p_cor),
         subtitle = paste0("P-value:", pv)) +
    xlab("TDH027 Average Gene #") +
    ylab("TDH027 S Score") +
    geom_label_repel(aes(label = Cluster),
                     size = 5,
                     col = "#3B3B53",
                     segment.color = "#3B3B53",
                     nudge_x = 0.03) +
    geom_smooth(method = lm, color="#AA4C26", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 40))
  ggsave(filename = paste0(outputDir, "TDH027_Correlation_Gene_Num_S_Score_NO_S12.png"), plot = p, width = 12, height = 10, dpi = 350)
  
  p_cor <- round(cor(plot_df$TDH027_Gene_Num,
                     plot_df$TDH027_G2M_Score, method = "spearman", use = "complete.obs"), 2)
  pv <- round(cor.test(plot_df$TDH027_Gene_Num,
                       plot_df$TDH027_G2M_Score, method = "spearman", use = "complete.obs")$p.value, 2)
  p <- ggplot(data = plot_df, aes(x=TDH027_Gene_Num, y=TDH027_G2M_Score)) +
    geom_point(col = "#487A8F", size = 8) +
    labs(title = paste0("Spearman Correlation:", p_cor),
         subtitle = paste0("P-value:", pv)) +
    xlab("TDH027 Average Gene #") +
    ylab("TDH027 G2M Score") +
    geom_label_repel(aes(label = Cluster),
                     size = 5,
                     col = "#3B3B53",
                     segment.color = "#3B3B53",
                     nudge_x = 0.03) +
    geom_smooth(method = lm, color="#AA4C26", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 40))
  ggsave(filename = paste0(outputDir, "TDH027_Correlation_Gene_Num_G2M_Score_NO_S12.png"), plot = p, width = 12, height = 10, dpi = 350)
  
  
  ### TDH019
  plot_df <- data.frame(Cluster=names(TDH019_clstr_gene_num)[-15],
                        TDH019_Gene_Num=TDH019_clstr_gene_num[-15],
                        TDH019_S_Score=TDH019_clstr_s_score[-15],
                        TDH019_G2M_Score=TDH019_clstr_g2m_score[-15],
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  p_cor <- round(cor(plot_df$TDH019_Gene_Num,
                     plot_df$TDH019_S_Score, method = "spearman", use = "complete.obs"), 2)
  pv <- round(cor.test(plot_df$TDH019_Gene_Num,
                       plot_df$TDH019_S_Score, method = "spearman", use = "complete.obs")$p.value, 2)
  p <- ggplot(data = plot_df, aes(x=TDH019_Gene_Num, y=TDH019_S_Score)) +
    geom_point(col = "#487A8F", size = 8) +
    labs(title = paste0("Spearman Correlation:", p_cor),
         subtitle = paste0("P-value:", pv)) +
    xlab("TDH019 Average Gene #") +
    ylab("TDH019 S Score") +
    geom_label_repel(aes(label = Cluster),
                     size = 5,
                     col = "#3B3B53",
                     segment.color = "#3B3B53",
                     nudge_x = 0.03) +
    geom_smooth(method = lm, color="#AA4C26", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 40))
  ggsave(filename = paste0(outputDir, "TDH019_Correlation_Gene_Num_S_Score_NO_S14.png"), plot = p, width = 12, height = 10, dpi = 350)
  
  p_cor <- round(cor(plot_df$TDH019_Gene_Num,
                     plot_df$TDH019_G2M_Score, method = "spearman", use = "complete.obs"), 2)
  pv <- round(cor.test(plot_df$TDH019_Gene_Num,
                       plot_df$TDH019_G2M_Score, method = "spearman", use = "complete.obs")$p.value, 2)
  p <- ggplot(data = plot_df, aes(x=TDH019_Gene_Num, y=TDH019_G2M_Score)) +
    geom_point(col = "#487A8F", size = 8) +
    labs(title = paste0("Spearman Correlation:", p_cor),
         subtitle = paste0("P-value:", pv)) +
    xlab("TDH019 Average Gene #") +
    ylab("TDH019 G2M Score") +
    geom_label_repel(aes(label = Cluster),
                     size = 5,
                     col = "#3B3B53",
                     segment.color = "#3B3B53",
                     nudge_x = 0.03) +
    geom_smooth(method = lm, color="#AA4C26", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 40))
  ggsave(filename = paste0(outputDir, "TDH019_Correlation_Gene_Num_G2M_Score_NO_S14.png"), plot = p, width = 12, height = 10, dpi = 350)
  
  
  
  ### load Jeremy's object
  JCC_Seurat_Obj <- readRDS(file = ExampleDataPath)
  
  ### see the UMAP plot
  DimPlot(object = JCC_Seurat_Obj, reduction = "umap",
          group.by = "seurat_clusters", label = TRUE,
          pt.size = 1.5)
  p <- FeaturePlot(JCC_Seurat_Obj, features = c("nFeature_RNA"),
                   cols = c("lightgray", "red"), label = TRUE, pt.size = 1.5, label.size = 10) +
    ggtitle("SJCAR19 Gene #") +
    theme_classic(base_size = 50) +
    theme(legend.key.size = unit(1, 'cm'))
  ggsave(paste0(outputDir, "SJCAR19_nFeature_UMAP.png"), plot = p, width = 15, height = 12, dpi = 350)
  
  
  
}

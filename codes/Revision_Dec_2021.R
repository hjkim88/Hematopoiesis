###
#   File name : Revision_Dec_2021.R
#   Author    : Hyunjin Kim
#   Date      : Dec 14, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : We received reviews from Nature Cell Biology. Revision should be done based on the comments.
#   
#   Instruction
#               1. Source("Revision_Dec_2021.R")
#               2. Run the function "trent_revision" - specify the input path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Revision_Dec_2021.R/Revision_Dec_2021.R")
#               > trent_revision(inputDataPath="./data/Combined_Seurat_Obj.RDS",
#                                outputDir="./results/Sep_2021/"
###

trent_revision <- function(inputDataPath="./data/Combined_Seurat_Obj.RDS",
                           outputDir="./results/Dec_2021/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  
  ## updated stroma RDS file
  Updated_Seurat_Obj <- readRDS(file = inputDataPath)
  
  ### create new column for the analysis
  Updated_Seurat_Obj@meta.data$Dev_Anno <- paste0(Updated_Seurat_Obj@meta.data$Development, "_",
                                                  Updated_Seurat_Obj@meta.data$Annotation)
  
  ### determine necessary variables
  time_points <- c("E16", "E18", "P0", "ADULT")
  HSPC_populations <- c("LTHSC", "STHSC", "MPP2", "MPP3", "MPP4")
  
  ### make the output directory
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  
  ### Reviewer #2 - 2.
  ### The annotation of clusters. The authors should provide clear information to readers how the clusters were assigned to the cell types.
  ### I didn’t see a single dotplot, violin plot, ridge plot, or heatmap that shows the level of expression of marker genes across all clusters.
  ### Instead authors use UMAPs which are not good for quantitative assessment of gene expression.
  ### In figure 4b it seems that authors drew by hand boundaries around various clusters which they divided by numbers.
  ### The relevance of these numbers remains unclear. For example, what is the difference between Mast1, Mast2i and Mast2ii, etc.?
  ### In supp. fig. 4b there is a number of UMAPs with expression of marker genes for each of the cell types.
  ### The expression bar is labelled low—high which has very little meaning in terms of quantification of expression.
  ### Furthermore, most of the marker genes are lowly expressed and often dispersed on the UMAP rather than confided to a distinct cluster
  ### e.g. s100a8, siglech etc.
  
  ### -> If you included the DE gene list of all the clusters, it must say those genes are differentially expressed in each of the cluster
  ### and the authors can look them up. But since the reviewer pointed this out, we can provide a dot plot of some marker genes
  ### for some important clusters.
  
  
  
  
  ### Reviewer #2 - 7.
  ### More detailed results for the scRNA-Seq analysis are required. The authors briefly describe the scRNA-Seq analysis in their
  ### Methods section, using PCA for dimensionality reduction, UMAP for visualisation and embedding, and Seurat FindAllMarkers
  ### to identify marker genes for each cluster. Further evidence of this work is required: 1) The Elbow curve to justify
  ### the choice of principal components; 2) Ridge/Violin plots of top marker genes for different clusters in the final set of annotations;
  ### 3) UMAP plots to indicate that there is no strong batch effect (no discussion of batch-correction software is provided in the Methods);
  ### 4) More details on the removal of doublets, including e.g. UMI distributions to indicate there are no bimodal effects,
  ### which might arise from doublets, in the data. As previously mentioned, dot plots in this study are misleading and should be
  ### rescaled based on something like the mean log(normalised gene expression).
  
  ### -> We can make the elbow plots, ridge/violin plots, UMAP plot of batch effect, UMI distribution etc.
  
  ### normalization
  Updated_Seurat_Obj <- NormalizeData(Updated_Seurat_Obj,
                                      normalization.method = "LogNormalize", scale.factor = 10000)
  ### find variable genes
  Updated_Seurat_Obj <- FindVariableFeatures(Updated_Seurat_Obj,
                                             selection.method = "vst", nfeatures = 2000)
  ### scaling
  Updated_Seurat_Obj <- ScaleData(Updated_Seurat_Obj,
                                  vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  ### PCA
  Updated_Seurat_Obj <- RunPCA(Updated_Seurat_Obj,
                               features = VariableFeatures(object = Updated_Seurat_Obj),
                               npcs = 50)
  ### UMAP
  Updated_Seurat_Obj <- RunUMAP(Updated_Seurat_Obj, dims = 1:50)
  
  ### 1. Elbow plot
  n <- 20
  
  ### PCA
  p <- ElbowPlot(Updated_Seurat_Obj, ndims = n, reduction = "pca") +
    ggtitle("Elbow Plot - PCA") +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30, color = "black", face = "bold"),
          axis.text = element_text(size = 25, color = "black", face = "bold"),
          axis.title = element_text(size = 30, color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'),
          legend.position = "right") +
    scale_x_continuous(labels = 1:n, breaks = 1:n)
  p$layers[[1]]$geom$default_aes$size <- 5
  ggsave(paste0(outputDir, "R2C7_Elbow_plot_PCA.png"), plot = p, width = 25, height = 10, dpi = 350)
  
  
  ### 3. batch effect
  ### Library, Time, Info, Tissue, Development, Cell_Type, HSPC
  
  ### variables
  vars <- c("Library", "Time", "Info", "Tissue", "Development", "Cell_Type", "HSPC", "Annotation")
  
  ### for each var, make a umap
  for(var in vars) {
    p <- DimPlot(object = Updated_Seurat_Obj, reduction = "umap", raster = TRUE,
                 group.by = var,
                 pt.size = 1) +
      ggtitle(paste0("UMAP - ", var)) +
      labs(color = var) +
      guides(colour = guide_legend(override.aes = list(size=10))) +
      theme_classic(base_size = 48) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 36, color = "black", face = "bold"),
            axis.title = element_text(size = 36, color = "black", face = "bold"),
            axis.text = element_text(size = 36, color = "black", face = "bold"),
            legend.title = element_text(size = 24, color = "black", face = "bold"),
            legend.text = element_text(size = 24, color = "black", face = "bold"),
            axis.ticks = element_blank())
    ggsave(paste0(outputDir, "UMAP_plot_", var, "_.png"), plot = p, width = 20, height = 10, dpi = 350)
    
    p <- DimPlot(object = Updated_Seurat_Obj, reduction = "pca", raster = TRUE,
                 group.by = var,
                 pt.size = 1) +
      ggtitle(paste0("PCA - ", var)) +
      labs(color = var) +
      guides(colour = guide_legend(override.aes = list(size=10))) +
      theme_classic(base_size = 48) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 36, color = "black", face = "bold"),
            axis.title = element_text(size = 36, color = "black", face = "bold"),
            axis.text = element_text(size = 36, color = "black", face = "bold"),
            legend.title = element_text(size = 24, color = "black", face = "bold"),
            legend.text = element_text(size = 24, color = "black", face = "bold"),
            axis.ticks = element_blank())
    ggsave(paste0(outputDir, "PCA_plot_", var, "_.png"), plot = p, width = 20, height = 10, dpi = 350)
  }
  
  
  ### 4. doublets; UMI distributions
  
  
  
  
  
  
  
  
  
}

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
  if(!require(ggthemes, quietly = TRUE)) {
    install.packages("ggthemes")
    require(ggthemes, quietly = TRUE)
  }
  remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
  require(DoubletFinder, quietly = TRUE)
  
  
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
  
  ### Cell cycle score (will be used later for regression out)
  Updated_Seurat_Obj <- CellCycleScoring(object = Updated_Seurat_Obj,
                                         g2m.features = cc.genes$g2m.genes,
                                         s.features = cc.genes$s.genes)
  
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
  
  
  ### Basic function to convert mouse to human gene names
  convertMouseGeneList <- function(x){
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    humanx <- unique(genesV2[, 2])
    # Print the first 6 genes found to the screen
    print(head(humanx))
    return(humanx)
  }
  
  ### Basic function to convert human to mouse gene names
  convertHumanGeneList <- function(x){
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
    humanx <- unique(genesV2[, 2])
    # Print the first 6 genes found to the screen
    print(head(humanx))
    return(humanx)
  }
  
  ### 4. doublets; UMI distributions
  ### https://github.com/chris-mcginnis-ucsf/DoubletFinder
  
  ### divide the data into different time points
  ### because I believe they were from different lanes
  E16_Seurat_Obj <- subset(Updated_Seurat_Obj,
                           cells = rownames(Updated_Seurat_Obj@meta.data)[which(Updated_Seurat_Obj$Development == "E16")])
  E18_Seurat_Obj <- subset(Updated_Seurat_Obj,
                           cells = rownames(Updated_Seurat_Obj@meta.data)[which(Updated_Seurat_Obj$Development == "E18")])
  P0_Seurat_Obj <- subset(Updated_Seurat_Obj,
                          cells = rownames(Updated_Seurat_Obj@meta.data)[which(Updated_Seurat_Obj$Development == "P0")])
  ADULT_Seurat_Obj <- subset(Updated_Seurat_Obj,
                             cells = rownames(Updated_Seurat_Obj@meta.data)[which(Updated_Seurat_Obj$Development == "ADULT")])
  
  ### preprocessing each seurat object
  seurat_preprocess <- function(so) {
    ### normalization
    so2 <- NormalizeData(so,
                         normalization.method = "LogNormalize", scale.factor = 10000)
    
    ### Cell cycle score (will be used later for regression out)
    so2 <- CellCycleScoring(object = so2,
                            g2m.features = convertHumanGeneList(cc.genes$g2m.genes),
                            s.features = convertHumanGeneList(cc.genes$s.genes))
    
    ### find variable genes
    so2 <- FindVariableFeatures(so2,
                                selection.method = "vst", nfeatures = 2000)
    ### scaling
    so2 <- ScaleData(so2,
                     vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
    ### PCA
    so2 <- RunPCA(so2,
                  features = VariableFeatures(object = so2),
                  npcs = 20)
    ### UMAP
    so2 <- RunUMAP(so2, dims = 1:20)
    
    return(so2)
  }
  
  E16_Seurat_Obj <- seurat_preprocess(E16_Seurat_Obj)
  E18_Seurat_Obj <- seurat_preprocess(E18_Seurat_Obj)
  P0_Seurat_Obj <- seurat_preprocess(P0_Seurat_Obj)
  ADULT_Seurat_Obj <- seurat_preprocess(ADULT_Seurat_Obj)
  
  ### function for doublet detection
    
  ### pK Identification (no ground-truth)
  sweep.res.list <- paramSweep_v3(so, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ### pk -> 0.01
  pk <- 0.01
  
  ### Homotypic Doublet Proportion Estimate
  ### Assuming 7.5% doublet formation rate
  nExp_poi <- round(0.075*nrow(so@meta.data))
  
  ### doublet detection with the parameters
  so <- doubletFinder_v3(so, PCs = 1:20, pK = pk, nExp = nExp_poi)
  
  
  
  
  ### Reviewer #2 - 2.
  ### The annotation of clusters. The authors should provide clear information to readers how the clusters were assigned to the cell types.
  ### I didn’t see a single dotplot, violin plot, ridge plot, or heatmap that shows the level of expression of marker genes across all clusters.
  ### Instead authors use UMAPs which are not good for quantitative assessment of gene expression. In figure 4b it seems that authors drew by hand
  ### boundaries around various clusters which they divided by numbers. The relevance of these numbers remains unclear. For example,
  ### what is the difference between Mast1, Mast2i and Mast2ii, etc.? In supp. fig. 4b there is a number of UMAPs with expression of marker genes
  ### for each of the cell types. The expression bar is labelled low—high which has very little meaning in terms of quantification of expression.
  ### Furthermore, most of the marker genes are lowly expressed and often dispersed on the UMAP rather than confided to a distinct cluster
  ### e.g. s100a8, siglech etc.
  
  ### -> Trent, I am going to make a dotplot to show gene expression changes among clusters. But could you provide few marker genes
  ### that you want to show that they are important for distinguishing some clusters? I would appreciate if you could also give me
  ### the order of the genes – the order you want to present in the plot.
  
  ### only use Heme data
  sub_seurat_obj <- subset(Updated_Seurat_Obj,
                           cells = rownames(Updated_Seurat_Obj@meta.data)[which(Updated_Seurat_Obj$Cell_Type == "Heme")])
  
  ### set genes & cell types of interest
  interesting_genes <- c("Bmi1", "Cd33", "Cd34", "Cd38", "Cd44", "Cd48", "Thy1", "Ly6a",
                         "Ly6e", "Myb", "Mcl1", "Pten", "Stat5a", "Stat5b")
  
  interesting_cell_types <- c("Mast1", "Mast2", "Mast3", "Mono1", "Prog", "Mono1", "Mono2", "Mono3", "TCell", "Bcell")
  
  ### only use interesting data
  sub_seurat_obj2 <- subset(sub_seurat_obj,
                            cells = rownames(sub_seurat_obj@meta.data)[which(sub_seurat_obj$Annotation %in% interesting_cell_types)])
  
  ### dot plot
  p <- DotPlot(sub_seurat_obj2,
               features = interesting_genes,
               cols = c("coral", "darkolivegreen"),
               group.by = "Annotation") +
    scale_size(range = c(2, 15)) +
    coord_flip() +
    xlab("") +
    ylab("") +
    theme_calc(base_size = 35) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file = paste0(outputDir, "R2C2_Dotplot_Markers.png"),
         plot = p, width = 20, height = 10, dpi = 350)
  
  ### ridge plot
  
  
  
  
  ### Reviewer #3 – 0.
  ### Could the early (E15.5-E18.5) FBM HSPCs are mainly CXCR4-? If yes, could the CXCR4- FBM HSPCs contribute to the increased
  ### lymphoid output (increased B-lymphoid output in transplantation assay E15.5-E18.5 (Fig 1c), increased T-lymphoid output in figure 6)?
  ### Could early BM niches are supportive for lymphoid biased HSCPs only, but not full-potential HSPCs?
    
  ### -> I can check if E16 & E18 are CXCR4- when compared to P0 & adult. We can discuss the further things after seeing the results.
  
  
  
  
  
  
  
}

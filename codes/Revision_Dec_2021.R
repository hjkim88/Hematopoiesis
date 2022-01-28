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
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  if(!require(data.table, quietly = TRUE)) {
    install.packages("data.table")
    require(data.table, quietly = TRUE)
  }
  if(!require(remotes, quietly = TRUE)) {
    install.packages("remotes")
    require(remotes, quietly = TRUE)
  }
  if(!require(reticulate, quietly = TRUE)) {
    install.packages("reticulate")
    require(reticulate, quietly = TRUE)
  }
  if(!require(RNAMagnet, quietly = TRUE)) {
    remotes::install_github("veltenlab/rnamagnet")
    require(RNAMagnet, quietly = TRUE)
  }
  if(!require(SeuratDisk, quietly = TRUE)) {
    remotes::install_github("mojaveazure/seurat-disk")
    require(SeuratDisk, quietly = TRUE)
  }
  # remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
  # require(DoubletFinder, quietly = TRUE)
  
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
  
  ### Elbow plot
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
  
  ### only use Heme data
  sub_seurat_obj <- subset(Updated_Seurat_Obj,
                           cells = rownames(Updated_Seurat_Obj@meta.data)[which(Updated_Seurat_Obj$Cell_Type == "Heme")])
  
  ### Elbow plot
  p <- ElbowPlot(sub_seurat_obj, ndims = n, reduction = "pca") +
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
  ggsave(paste0(outputDir, "R2C7_Elbow_plot_PCA2.png"), plot = p, width = 25, height = 10, dpi = 350)
  
  
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
  
  
  
  
  ### Reviewer #1 - 4.
  ### Beaudin et al. defined a Flt3+ HSC population that is transient. Is this population evident in the authors dataset?
  
  ### -> I think I need to talk with Trent to know the meaning of "transient" of the Flt3+ HSC population,
  ### but this is definitely something that we can do with our dataset.
  
  ### UMAP of Flt3+ population
  p <- FeaturePlot(Updated_Seurat_Obj, features = "Flt3",
                   cols = c("darkolivegreen", "coral"),
                   raster = TRUE,
                   pt.size = 1,
                   label = TRUE,
                   label.size = 4,
                   ncol = 1)
  
  p2 <- DimPlot(object = Updated_Seurat_Obj, reduction = "umap", raster = TRUE,
                group.by = "Development",
                pt.size = 1)
  
  p3 <- DimPlot(object = Updated_Seurat_Obj, reduction = "umap", raster = TRUE,
                group.by = "HSPC",
                pt.size = 1)
  
  g <- arrangeGrob(grobs = list(p, p2, p3),
                   nrow = 2,
                   ncol = 2,
                   top = "")
  ggsave(paste0(outputDir, "R1C4_UMAP_Flt3.png"), plot = g, width = 12, height = 10, dpi = 350)
  
  
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
  p <- RidgePlot(sub_seurat_obj2,
                 features = interesting_genes,
                 group.by = "Annotation",
                 ncol = 5)
  for(i in 1:length(interesting_genes)) {
    p[[i]] <- p[[i]] +
      ylab("") +
      theme_classic(base_size = 10) +
      theme(legend.key.size = unit(3, 'cm'),
            legend.position = "none",
            legend.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
            legend.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35, color = "black", face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, size = 25, color = "black", face = "bold"),
            axis.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
            axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"))
  }
  
  ggsave(file = paste0(outputDir, "R2C2_Ridgeplot_Markers.png"),
         plot = p, width = 25, height = 15, dpi = 350)
  
  
  ### violin plot
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2$Annotation)
  p <- VlnPlot(sub_seurat_obj2, features = interesting_genes, slot = "data",
               pt.size = 0, ncol = 4)
  for(i in 1:length(interesting_genes)) {
    p[[i]] <- p[[i]] + geom_boxplot(width=0.1) +
      stat_compare_means(size = 5) +
      ylab("Log-Normalized EXP") +
      stat_summary(fun=mean, geom="point", size=3, color="#852A30") +
      theme_classic(base_size = 15) +
      theme(legend.position = "none",
            axis.title.x = element_blank())
  }
  
  ### save the violin plot
  ggsave(file = paste0(outputDir, "R2C2_Violin_Markers.png"),
         plot = p, width = 30, height = 20, dpi = 350)
  
  
  ### Reviewer #3 – 0.
  ### Could the early (E15.5-E18.5) FBM HSPCs are mainly CXCR4-? If yes, could the CXCR4- FBM HSPCs contribute to the increased
  ### lymphoid output (increased B-lymphoid output in transplantation assay E15.5-E18.5 (Fig 1c), increased T-lymphoid output in figure 6)?
  ### Could early BM niches are supportive for lymphoid biased HSCPs only, but not full-potential HSPCs?
  
  ### -> I can check if E16 & E18 are CXCR4- when compared to P0 & adult. We can discuss the further things after seeing the results.
  
  ### dot plot
  p <- DotPlot(sub_seurat_obj,
               features = "Cxcr4",
               cols = c("coral", "darkolivegreen"),
               group.by = "Development") +
    scale_size(range = c(10, 30)) +
    coord_flip() +
    xlab("") +
    ylab("") +
    theme_classic(base_size = 35) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file = paste0(outputDir, "R3C0_Dotplot_Cxcr4.png"),
         plot = p, width = 18, height = 12, dpi = 350)
  
  ### violin plot
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj$Development)
  p <- VlnPlot(sub_seurat_obj, features = "Cxcr4", slot = "data",
               pt.size = 0, ncol = 1)
  
  p <- p + geom_boxplot(width=0.1) +
    stat_compare_means(size = 5) +
    ylab("Log-Normalized EXP") +
    stat_summary(fun=mean, geom="point", size=3, color="#852A30") +
    theme_classic(base_size = 35) +
    theme(legend.position = "none",
          axis.title.x = element_blank())
  
  ### save the violin plot
  ggsave(file = paste0(outputDir, "R3C0_Violin_Cxcr4.png"),
         plot = p, width = 15, height = 10, dpi = 350)
  
  
  ### Reviewer #3 – 5.
  ### In figure 5 a-c, the author annotated 3 different E16.5 FBM LT-HSCs showing different biological process with GO term analysis.
  ### The author should go further to analysis the expression profile of HSCs signature genes, cell cycle related genes and chemotaxis
  ### related genes, e.g., SLAMF1, PROCR, HOXA9, CD48, CXCR4, CXCR3, ACKR7 receptors. Is it possible that not all of the 3 different clusters
  ### are full functional, unbiased LT-HSCs?
  
  ### -> We might want to show a dot plot of the genes that the reviewer mentioned to compare the expression profiles among the 3 clusters.
  ### Then we can additionally describe what are the differences of the 3 clusters.
  
  ### it's actually about LTHSCs - 3 clusters of E16.5
  LTHSC_E16_Seurat_Obj <- subset(E16_Seurat_Obj,
                                 cells = rownames(E16_Seurat_Obj@meta.data)[which(E16_Seurat_Obj$HSPC == "LTHSC")])
  
  ### draw UMAP to find the 3 clusters
  
  
  
  
  
  
  ### Reviewer #1 - 3.
  ### The authors should also incorporate some comparative analyses with the human fetal bone marrow recently published
  ### by Jardine et al. It would be helpful for the field to have some high level discussion of similarities and differences.
  
  ### -> 1.) Pull out what they have annotated as stroma and HSC/MPP and run them through our RNA-Magnet pipeline
  ### (with ligand coming from Stroma and receptor coming from HSC/MPP).
  ### 2.) Pull out what they have annotated as HSC/MPP, recluster, and then perform GO analysis on the clusters
  ### to see if there is any overlap in the GO terms between their HSC/MPPs and our HSCs.
  
  # Basic function to convert mouse to human gene names
  convertMouseGeneList <- function(x){
    if(!require(biomaRt, quietly = TRUE)) {
      BiocManager::install("biomaRt")
      require(biomaRt, quietly = TRUE)
    }
    
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    
    return(genesV2)
  }
  
  # Basic function to convert human to mouse gene names
  convertHumanGeneList <- function(x){
    if(!require(biomaRt, quietly = TRUE)) {
      BiocManager::install("biomaRt")
      require(biomaRt, quietly = TRUE)
    }
    
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
    
    return(genesV2)
  }
  
  #'Low-level function to run core RNA magnet steps
  #'
  #' Users are advised to use the top-level functions \code{\link{RNAMagnetAnchors}} and \code{\link{RNAMagnetSignaling}} which appropriately set default parameters and return user-friendly return values.
  #' This is a low level function for development purposes.
  #'@param seurat An object of class \code{\link[Seurat]{seurat}} containing a valid clustering and t-SNE information
  #'
  #'. For information on how to create such an object, see https://satijalab.org/seurat/get_started.html
  #'@param anchors A character vector of anchor populations. Entries must be levels of the seurat identity vector. If \code{NULL}: All entries of the seurat identity vector are used as anchors.
  #'@param neighborhood.distance See \code{\link{RNAMagnetAnchors}}
  #'@param neighborhood.gradient See \code{\link{RNAMagnetAnchors}}
  #'@param .k Fuzzification parameter, see detail. Recommended to leave at the default value.
  #'@param .x0 Fuzzification parameter, see detail. Recommended to leave at the default value.
  #'@param .minExpression Minimal expression level of genes to be included, specified as the number of cells in the dataset that express the gene.
  #'@param .minMolecules Number of molecules per cell required to count the cells as positive.
  #'@param .version The version of the underlying ligand-receptor database. See \code{\link{getLigandsReceptors}}.
  #'@param .cellularCompartment Types of ligands to be included. For physical interactions, defaults to \code{ c("Membrane","ECM","Both")}.  See \code{\link{getLigandsReceptors}}.
  #'@param .manualAnnotation Annotation status of ligands to be included. Default to \code{"Correct"}. See \code{\link{getLigandsReceptors}}.
  #'@param .symmetric Assume that if A is a receptor for B, B is also a receptor for A
  #'@details The algorithm takes the following steps: \enumerate{
  #'\item Ligand-receptor pairs are selected based on the parameters \code{.version}, \code{.cellularCompartment} and \code{.manualAnnotation}. Choice of \code{.cellularCompartment} is crucial for determining the algorithm's behavior, e.g. if set to \code{c("Secreted","Both")}, paracrine signaling interactions involving soluble ligands are investigated.
  #'\item Dropout values in the expression levels of ligands and receptors are imputed using \code{\link[Rmagic]{magic}}
  #'\item Mean expression level of ligands and receptors is computed for all anchor populations
  #'\item For each cell or anchor population, the expression of each ligand and receptor is encoded as a fuzzy logic variable
  #'\item Fuzzy logic AND is used to compute probabilities for a given interaction to be active between a single cell and an anchor population
  #'\item An interaction score is computed as the sum of interaction probabilities across all possible ligand-receptor pairs
  #'\item Specificty scores are computed by comparing interaction scores to average interaction scores in a local neighborhood.
  #'}
  #'@details Add the methods section of the paper here!
  #'@return Returns an object of class \code{\link{rnamagnet}}
  #'@export
  RNAMagnetBase <- function(seurat, anchors=NULL, human=FALSE, neighborhood.distance=NULL, neighborhood.gradient =NULL, .k = 10, .x0 = 0.5, .minExpression,.minMolecules=1, .version = "1.0.0", .cellularCompartment, .manualAnnotation = "Correct", .symmetric = F) {
    cat("Setting everything up...\n")
    
    if (grepl("^3", Biobase::package.version("Seurat")) || grepl("^4", Biobase::package.version("Seurat"))) {
      seurat.ident <- Idents(seurat)
      seurat.cell.names <- colnames(seurat)
      seurat.pca <- Embeddings(seurat, reduction = "pca")
      seurat.raw.data <- GetAssayData(seurat, slot = "counts")
    } else {
      seurat.ident <- seurat@ident
      seurat.cell.names <- seurat@cell.names
      seurat.pca <- seurat@dr$pca@cell.embeddings
      seurat.raw.data <- seurat@raw.data
    }
    
    ### if the input seurat is a human data
    if(human) {
      ### get mapping data
      tempGenes <- convertHumanGeneList(rownames(seurat.raw.data))
      
      ### only retain genes those have mouse genes
      seurat.raw.data <- seurat.raw.data[unique(tempGenes[,1]),]
      
      ### get human gene duplicates
      dups <- unique(tempGenes[,1][which(duplicated(tempGenes[,1]))])
      nonduplicated_tempGenes <- tempGenes[which(!(tempGenes[,1] %in% dups)),]
      duplicated_tempGenes <- tempGenes[which(tempGenes[,1] %in% dups),]
      
      ### only retain genes those have mouse genes
      seurat.raw.data2 <- seurat.raw.data[nonduplicated_tempGenes[,1],]
      rownames(seurat.raw.data2) <- nonduplicated_tempGenes[,2]
      seurat.raw.data3 <- seurat.raw.data[unique(duplicated_tempGenes[,1]),]
      
      ### one human gene can have multiple mouse genes
      ### copy all the gexp for those multiple mouse genes
      for(i in 1:nrow(seurat.raw.data3)) {
        temp <- tempGenes[which(tempGenes[,1] == rownames(seurat.raw.data3)[i]),]
        rownames(seurat.raw.data3)[i] <- temp[1,2]
        for(j in 1:(nrow(temp)-1)) {
          seurat.raw.data3 <- rbind(seurat.raw.data3, seurat.raw.data3[i,])
        }
        rownames(seurat.raw.data3)[(nrow(seurat.raw.data3)-nrow(temp)+2):nrow(seurat.raw.data3)] <- temp[2:nrow(temp),2]
      }
      
      ### combine the partial results
      seurat.raw.data <- rbind(seurat.raw.data2, seurat.raw.data3)
      
      ### check if adding multiple mouse genes were done accurately
      print(identical(nrow(seurat.raw.data), nrow(tempGenes)))
      print(identical(length(which(duplicated(rownames(seurat.raw.data)))), length(which(duplicated(tempGenes[,2])))))
      
      ### now there are duplicated mouse genes
      ### get mean value for those
      dup_idx <- which(duplicated(rownames(seurat.raw.data)))
      dups <- unique(rownames(seurat.raw.data)[dup_idx])
      for(dup in dups) {
        target_idx <- which(rownames(seurat.raw.data) == dup)
        seurat.raw.data[target_idx[1],] <- apply(seurat.raw.data[target_idx,], 2, median)
      }
      seurat.raw.data <- seurat.raw.data[-dup_idx,]
    }
    
    if (is.null(anchors)) anchors <- as.character(unique(seurat.ident))
    
    out <- new("rnamagnet", celltype = seurat.ident, params = list("neighborhood.distance"=neighborhood.distance, "neighborhood.gradient" =neighborhood.gradient, ".k" = .k, ".x0" = .x0, ".minExpression" = .minExpression, ".minMolecules" = .minMolecules, ".cellularCompartment" = .cellularCompartment, ".manualAnnotation" = .manualAnnotation, ".symmetric" = .symmetric))
    
    #compute cell-cell similarity
    similarity <- as.matrix(1-cor(t(seurat.pca[,1:15])))
    
    #prepare database
    ligrec <- getLigandsReceptors(.version, .cellularCompartment, .manualAnnotation)
    if (.symmetric) ligrec <- makeSymmetric(ligrec) #for physical interactions: i a binds b, b binds a.
    
    #prepare genes included into MAGIC
    filteredGenes <- rownames(seurat.raw.data)[apply(seurat.raw.data[,seurat.cell.names] >=.minMolecules ,1,sum) > .minExpression]
    
    genes <- unique(c(ligrec$Receptor.Mouse,
                      ligrec$Ligand.Mouse))
    
    genes <- genes[sapply(genes, function(x) {
      entries <- strsplit(x, "[&|]")
      if (grepl("&", x))
        all(entries[[1]] %in% filteredGenes)
      else any(entries[[1]] %in% filteredGenes)
    })]
    
    genes_formagic <- unlist(strsplit( genes,"[&|]"))
    
    
    formagic <- Matrix::t(seurat.raw.data[,seurat.cell.names])
    formagic <- Rmagic::library.size.normalize(formagic)
    formagic <- sqrt(formagic)
    
    #Run MAGIC
    cat("Now running MAGIC to impute dropout values...\n")
    mymagic <- Rmagic::magic(formagic, genes = genes_formagic, seed = 0xbeef)
    mymagic <- as.data.frame(mymagic)
    
    #handle expression levels of heterodimers
    resolvedRawData <- resolveData(t(mymagic), genes)
    resolvedRawData <- t(apply(resolvedRawData, 1, function(x) (x - min(x )) / (max(x) - min(x)) ))
    
    
    annotated_genes <- rownames(resolvedRawData)
    out@mylr <- subset(ligrec, Receptor.Mouse %in% annotated_genes &
                         Ligand.Mouse %in% annotated_genes)
    
    stepf<- Vectorize(function(x) if (x<0) 0 else x)
    
    cat("Now running RNAMagnet...\n")
    
    #compute mean gene expression level per population
    out@anchors <- do.call(cbind, lapply(anchors, function(id) {
      apply(resolvedRawData[,seurat.ident == id,drop=FALSE],1,mean)
    }));
    colnames(out@anchors) <- anchors
    
    #use fuzzy logic and operations to compute interaction score
    out@interaction <- sapply(anchors, function(pop_l) {
      out@mylr$expression_ligand <- out@anchors[out@mylr$Ligand.Mouse, pop_l]
      sapply(seurat.cell.names, function(cell_r) {
        out@mylr$expression_receptor <-resolvedRawData[out@mylr$Receptor.Mouse, cell_r]
        sum(kernel(out@mylr$expression_ligand, k =.k, x0 = .x0) * kernel(out@mylr$expression_receptor, k = .k, x0=.x0 )) #kernel performs fuzzification
      })
    })
    ### sometimes this is missing
    colnames(out@interaction) <- anchors
    
    out@specificity <- t(sapply(rownames(out@interaction), function(cell) {
      x <- out@interaction[cell,]
      beta <- x/sum(x)
      if (!is.null(neighborhood.distance)) alpha <- apply(out@interaction * (1-kernel(similarity[cell,],neighborhood.gradient,x0=neighborhood.distance )),2,sum) else alpha <- apply(out@interaction,2,mean)
      alpha <- alpha / sum(alpha)
      beta- alpha
    }))
    rownames(out@specificity) <- rownames(out@interaction)
    
    out@adhesiveness <- apply(mymagic, 1,function(x) sum(kernel(x,x0 = .x0, k = .k)))
    return(out)
    
  }
  
  ### RNAMagnetAnchors
  RNAMagnetAnchors <- function(seurat, anchors, human=FALSE, return = "summary", neighborhood.distance = 0.7, neighborhood.gradient = 3, .k = 10, .x0 = 0.5, .minExpression = 0, .minMolecules = 1, .version = "1.0.0", .cellularCompartment = c("Membrane","ECM","Both"), .manualAnnotation = "Correct" ) {
    
    myMagnet <- RNAMagnetBase(seurat, anchors, human, neighborhood.distance,neighborhood.gradient, .k, .x0, .minExpression, .minMolecules, .version, .cellularCompartment, .manualAnnotation,TRUE)
    if (return=="rnamagnet-class") myMagnet else data.frame(direction = as.factor(colnames(myMagnet@specificity)[apply(myMagnet@specificity,1,which.max)]), adhesiveness = myMagnet@adhesiveness, myMagnet@specificity[,anchors])
    
  }
  
  ### RNAMagnetSignaling
  RNAMagnetSignaling <- function(seurat, human=FALSE, neighborhood.distance = NULL, neighborhood.gradient = NULL, .k = 10, .x0 = 0.5, .minExpression = 10, .minMolecules = 1, .version = "1.0.0", .cellularCompartment = c("Secreted","Both"), .manualAnnotation = "Correct" ) {
    
    RNAMagnetBase(seurat, anchors = NULL, human, neighborhood.distance,neighborhood.gradient, .k, .x0, .minExpression, .minMolecules, .version, .cellularCompartment, .manualAnnotation, FALSE)
    
  }
  
  # ### load the Jardine data
  # jardine_meta <- read.table(file = "./data/E-MTAB-9389.sdrf.txt",
  #                            sep = "\t", header = TRUE,
  #                            stringsAsFactors = FALSE, check.names = FALSE)
  # jardine_data <- fread(file = "./data/fbm_tot_gex_20201104.csv",
  #                       header = TRUE, 
  #                       stringsAsFactors = FALSE, check.names = FALSE)
  # 
  # ### change data.table to data.frame
  # jardine_data <- setDF(jardine_data)
  # 
  # ### set gene names as row names
  # rownames(jardine_data) <- jardine_data[,1]
  # jardine_data <- jardine_data[,-1]
  # 
  # ### make a seurat object from the count matrix
  # jardine_seurat <- CreateSeuratObject(counts = jardine_data,
  #                                      project = "Jardine",
  #                                      assay = "RNA",
  #                                      meta.data = NULL)
  # 
  # ### pre-process
  # 
  # ### MT percentage
  # jardine_seurat[["percent.mt"]] <- PercentageFeatureSet(jardine_seurat, pattern = "^MT-")
  # 
  # ### normalization
  # jardine_seurat <- NormalizeData(jardine_seurat,
  #                                 normalization.method = "LogNormalize", scale.factor = 10000)
  # 
  # ### Cell cycle score (will be used later for regression out)
  # jardine_seurat <- CellCycleScoring(object = jardine_seurat,
  #                                    g2m.features = cc.genes$g2m.genes,
  #                                    s.features = cc.genes$s.genes)
  # 
  # ### find variable genes
  # jardine_seurat <- FindVariableFeatures(jardine_seurat,
  #                                        selection.method = "vst", nfeatures = 2000)
  # 
  # ### scaling
  # jardine_seurat <- ScaleData(jardine_seurat,
  #                             vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  # 
  # ### PCA
  # jardine_seurat <- RunPCA(jardine_seurat,
  #                          features = VariableFeatures(object = jardine_seurat),
  #                          npcs = 50)
  # 
  # ### draw elbow plot to see how many PCAs are appropriate
  # ElbowPlot(jardine_seurat, ndims = 50, reduction = "pca")
  # 
  # ### UMAP
  # jardine_seurat <- RunUMAP(jardine_seurat, dims = 1:50)
  # 
  # ### save the object
  # saveRDS(object = jardine_seurat,
  #         file = "./data/Jardine_SeuratObj.RDS")
  # 
  # ### print UMAP
  # DimPlot(object = jardine_seurat, reduction = "umap", raster = TRUE,
  #         group.by = NULL, pt.size = 1)
  
  ### we got a combined and annotated data so load it instead of building from scratch
  
  ### convert h5ad to h5seurat and load
  # Convert(source = "./data/fig1b_fbm_scaled_gex_updated_dr_20210104.h5ad", dest = "h5seurat")
  jardine_seurat2 <- LoadH5Seurat("./data/fig1b_fbm_scaled_gex_updated_dr_20210104.h5seurat")
  
  # ### https://satijalab.org/seurat/archive/v2.4/conversion_vignette.html
  # ad <- import("anndata", convert = FALSE)
  # obj <- ad$read_h5ad("./data/fig1b_fbm_scaled_gex_updated_dr_20210104.h5ad")
  # jardine_seurat3 <- Convert(obj, to = "seurat")
  
  ### see UMAP of the h5ad file
  DimPlot(object = jardine_seurat2, reduction = "umap", raster = TRUE,
          group.by = "broad_fig1_cell.labels", pt.size = 1)
  
  ### pull out HSC/MPP and Stroma only
  print(levels(jardine_seurat2$broad_fig1_cell.labels))
  subset_jardine_obj <- subset(jardine_seurat2,
                               cells = rownames(jardine_seurat2@meta.data)[which(jardine_seurat2$broad_fig1_cell.labels %in% c("HSC/MPP and pro", "stroma"))])
  
  ### PCA is not found so we run PCA
  subset_jardine_obj <- FindVariableFeatures(subset_jardine_obj)
  subset_jardine_obj <- ScaleData(subset_jardine_obj)
  subset_jardine_obj <- RunPCA(subset_jardine_obj, npcs = 15)
  
  ### run RNAMagnet with anchors
  subset_jardine_obj <- SetIdent(object = subset_jardine_obj,
                                 cells = rownames(subset_jardine_obj@meta.data),
                                 value = subset_jardine_obj$broad_fig1_cell.labels)
  result <- RNAMagnetAnchors(seurat = subset_jardine_obj,
                             anchors = unique(subset_jardine_obj$broad_fig1_cell.labels),
                             human = TRUE)
  
  
  
  
  
  
}

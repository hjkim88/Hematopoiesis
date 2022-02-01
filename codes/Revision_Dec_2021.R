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
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  if(!require(org.Mm.eg.db, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Mm.eg.db")
    require(org.Mm.eg.db, quietly = TRUE)
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
  
  ### UMI distribution
  temp <- apply(Updated_Seurat_Obj@assays$RNA@counts, 2, mean)
  
  
  
  ### 4. doublets; UMI distributions
  ### https://github.com/chris-mcginnis-ucsf/DoubletFinder
  
  ### Basic function to convert mouse to human gene names
  convertMouseGeneList <- function(x){
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    humanx <- unique(genesV2[, 2])
    
    return(humanx)
  }
  
  ### Basic function to convert human to mouse gene names
  convertHumanGeneList <- function(x){
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
    humanx <- unique(genesV2[, 2])
    
    return(humanx)
  }
  
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
  ### Jardine: 13-17 week human data
  
  ### -> 1.) Pull out what they have annotated as stroma and HSC/MPP and run them through our RNA-Magnet pipeline
  ### (with ligand coming from Stroma and receptor coming from HSC/MPP).
  ### 2.) Pull out what they have annotated as HSC/MPP, recluster, and then perform GO analysis on the clusters
  ### to see if there is any overlap in the GO terms between their HSC/MPPs and our HSCs.
  
  # Basic function to convert mouse to human gene names
  convertMouseGeneList2 <- function(x){
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
  convertHumanGeneList2 <- function(x){
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
      tempGenes <- convertHumanGeneList2(rownames(seurat.raw.data))
      
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
  RNAMagnetSignaling <- function(seurat, anchors=NULL, human=FALSE, neighborhood.distance = NULL, neighborhood.gradient = NULL, .k = 10, .x0 = 0.5, .minExpression = 10, .minMolecules = 1, .version = "1.0.0", .cellularCompartment = c("Secreted","Both"), .manualAnnotation = "Correct" ) {
    
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
  ElbowPlot(subset_jardine_obj, ndims = 15, reduction = "pca")
  subset_jardine_obj <- RunUMAP(subset_jardine_obj, dims = 1:15)
  
  ### run RNAMagnet with anchors
  subset_jardine_obj <- SetIdent(object = subset_jardine_obj,
                                 cells = rownames(subset_jardine_obj@meta.data),
                                 value = subset_jardine_obj$broad_fig1_cell.labels)
  subset_jardine_obj$broad_fig1_cell.labels <- factor(subset_jardine_obj$broad_fig1_cell.labels,
                                                      levels = as.character(unique(subset_jardine_obj$broad_fig1_cell.labels)))
  result <- RNAMagnetAnchors(seurat = subset_jardine_obj,
                             anchors = levels(subset_jardine_obj$broad_fig1_cell.labels),
                             human = TRUE)
  # saveRDS(result, file = "./rnamagnet_anchor_result.rds")
  
  ### write the result as an Excel file
  write.xlsx2(data.frame(Cell=rownames(result), result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "RNAMagnet_Anchors_Jardine.xlsx"),
              sheetName = paste0("RNAMagnet_Result"),
              row.names = FALSE)
  
  ### add RNAMagnet info to the seurat object
  subset_jardine_obj$direction <- as.character(result[rownames(subset_jardine_obj@meta.data),"direction"])
  subset_jardine_obj$adhesiveness <- as.numeric(result[rownames(subset_jardine_obj@meta.data),"adhesiveness"])
  subset_jardine_obj$specificity <- as.numeric(sapply(rownames(subset_jardine_obj@meta.data), function(x) {
    result[x, result[x,"direction"]]
  }))
  
  ### make a data frame for ggplot
  dim_map <- Embeddings(subset_jardine_obj, reduction = "umap")[rownames(subset_jardine_obj@meta.data),]
  plot_df <- data.frame(X=dim_map[rownames(subset_jardine_obj@meta.data),1],
                        Y=dim_map[rownames(subset_jardine_obj@meta.data),2],
                        cluster_color = subset_jardine_obj$broad_fig1_cell.labels,
                        direction = subset_jardine_obj@meta.data$direction,
                        adhesiveness = subset_jardine_obj@meta.data$adhesiveness,
                        specificity = subset_jardine_obj@meta.data$specificity,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### draw a scatter plot with the adhesiveness info
  p <- list()
  p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
    xlab("UMAP1") + ylab("UMAP2") +
    labs(col="Cluster") +
    ggtitle("UMAP with Cell Type") +
    scale_color_brewer(palette="Dark2") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'))
  ### the id column in the plot_df should be a factor
  # p[[1]] <- LabelClusters(plot = p[[1]], id = "cluster_color", col = "black")
  
  p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="direction", alpha="adhesiveness"), size=2) +
    xlab("UMAP1") + ylab("UMAP2") +
    labs(col="Direction", alpha="Adhesiveness") +
    ggtitle("UMAP with Direction & Adhesiveness") +
    scale_color_brewer(palette="Set1") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'))
  ### the id column in the plot_df should be a factor
  # p[[2]] <- LabelClusters(plot = p[[2]], id = "cluster_color", col = "black")
  
  p[[3]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="direction", alpha="specificity"), size=2) +
    xlab("UMAP1") + ylab("UMAP2") +
    labs(col="Direction", alpha="Specificity") +
    ggtitle("UMAP with Direction & Specificity") +
    scale_color_brewer(palette="Set1") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'))
  ### the id column in the plot_df should be a factor
  # p[[3]] <- LabelClusters(plot = p[[3]], id = "cluster_color", col = "black")
  
  ### save the plots
  g <- arrangeGrob(grobs = list(p[[1]], p[[3]]),
                   nrow = 1,
                   ncol = 2,
                   top = "")
  ggsave(file = paste0(outputDir, "RNAMagnet_Anchors_Jardine.png"), g, width = 30, height = 10, dpi = 350)
  
  
  ### RNAMagnet signaling
  result2 <- RNAMagnetSignaling(seurat = subset_jardine_obj,
                                human = TRUE)
  # saveRDS(result2, file = "./rnamagnet_signaling_result.rds")
  
  ### get all the interaction list
  interaction_list <- NULL
  for(clust1 in levels(subset_jardine_obj$broad_fig1_cell.labels)) {
    for(clust2 in levels(subset_jardine_obj$broad_fig1_cell.labels)) {
      il <- getRNAMagnetGenes(result2, clust1, clust2, thresh = 0)
      if(nrow(il) > 0) {
        il$ligand_cluster <- clust1
        il$receptor_cluster <- clust2
        if(is.null(interaction_list)) {
          interaction_list <- il
        } else {
          interaction_list <- rbind(interaction_list, il)
        }
      }
    }
  }
  interaction_list <- cbind(interaction_list[3:4], interaction_list[1:2])
  colnames(interaction_list) <- c("Ligand_Cluster", "Receptor_Cluster", "Interaction_Score", "Interaction_Pair")
  
  ### write out the interaction list
  write.xlsx2(interaction_list,
              file = paste0(outputDir, "RNAMagnet_Signaling_Jardine.xlsx"),
              sheetName = paste0("Signaling_RNAMagnet_Interaction_List"),
              row.names = FALSE)
  
  ###
  ### Now pull ours
  ###
  
  ### define our annotation specifically
  ### MPP/Stroma & E16, E18, and P0
  Updated_Seurat_Obj$New_HSPC <- Updated_Seurat_Obj$HSPC
  Updated_Seurat_Obj$New_HSPC[which(Updated_Seurat_Obj$New_HSPC %in% c("LTHSC", "STHSC", "MPP2", "MPP3", "MPP4"))] <- "HSPC"
  Updated_Seurat_Obj$New_Anno <- paste0(Updated_Seurat_Obj$Development, "_", Updated_Seurat_Obj$New_HSPC)
  
  ### pull out E16 HSC/MPP and Stroma only
  print(unique(Updated_Seurat_Obj$New_Anno))
  subset_our_obj1 <- subset(Updated_Seurat_Obj,
                            cells = rownames(Updated_Seurat_Obj@meta.data)[which(Updated_Seurat_Obj$New_Anno %in% c("E16_HSPC", "E16_Stroma"))])
  
  ### run PCA and UMAP
  subset_our_obj1 <- FindVariableFeatures(subset_our_obj1)
  subset_our_obj1 <- ScaleData(subset_our_obj1)
  subset_our_obj1 <- RunPCA(subset_our_obj1, npcs = 15)
  ElbowPlot(subset_our_obj1, ndims = 15, reduction = "pca")
  subset_our_obj1 <- RunUMAP(subset_our_obj1, dims = 1:15)
  
  ### run RNAMagnet with anchors
  subset_our_obj1 <- SetIdent(object = subset_our_obj1,
                              cells = rownames(subset_our_obj1@meta.data),
                              value = subset_our_obj1$New_Anno)
  subset_our_obj1$New_Anno <- factor(subset_our_obj1$New_Anno,
                                     levels = as.character(unique(subset_our_obj1$New_Anno)))
  our_result_E16 <- RNAMagnetAnchors(seurat = subset_our_obj1,
                                     anchors = levels(subset_our_obj1$New_Anno),
                                     human = FALSE)
  
  ### write the result as an Excel file
  write.xlsx2(data.frame(Cell=rownames(our_result_E16), our_result_E16,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "RNAMagnet_Anchors_Ours_E16.xlsx"),
              sheetName = paste0("RNAMagnet_Result"),
              row.names = FALSE)
  
  ### add RNAMagnet info to the seurat object
  subset_our_obj1$direction <- as.character(our_result_E16[rownames(subset_our_obj1@meta.data),"direction"])
  subset_our_obj1$adhesiveness <- as.numeric(our_result_E16[rownames(subset_our_obj1@meta.data),"adhesiveness"])
  subset_our_obj1$specificity <- as.numeric(sapply(rownames(subset_our_obj1@meta.data), function(x) {
    our_result_E16[x, our_result_E16[x,"direction"]]
  }))
  
  ### make a data frame for ggplot
  dim_map <- Embeddings(subset_our_obj1, reduction = "umap")[rownames(subset_our_obj1@meta.data),]
  plot_df <- data.frame(X=dim_map[rownames(subset_our_obj1@meta.data),1],
                        Y=dim_map[rownames(subset_our_obj1@meta.data),2],
                        cluster_color = subset_our_obj1$New_Anno,
                        direction = subset_our_obj1@meta.data$direction,
                        adhesiveness = subset_our_obj1@meta.data$adhesiveness,
                        specificity = subset_our_obj1@meta.data$specificity,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### draw a scatter plot with the adhesiveness info
  p <- list()
  p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
    xlab("UMAP1") + ylab("UMAP2") +
    labs(col="Cluster") +
    ggtitle("UMAP with Cell Type") +
    scale_color_brewer(palette="Dark2") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'))
  ### the id column in the plot_df should be a factor
  # p[[1]] <- LabelClusters(plot = p[[1]], id = "cluster_color", col = "black")
  
  p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="direction", alpha="specificity"), size=2) +
    xlab("UMAP1") + ylab("UMAP2") +
    labs(col="Direction", alpha="Specificity") +
    ggtitle("UMAP with Direction & Specificity") +
    scale_color_brewer(palette="Set1") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'))
  ### the id column in the plot_df should be a factor
  # p[[2]] <- LabelClusters(plot = p[[3]], id = "cluster_color", col = "black")
  
  ### save the plots
  g <- arrangeGrob(grobs = p,
                   nrow = 1,
                   ncol = 2,
                   top = "")
  ggsave(file = paste0(outputDir, "RNAMagnet_Anchors_Ours_E16.png"), g, width = 30, height = 10, dpi = 350)
  
  
  ### pull out E18 HSC/MPP and Stroma only
  print(unique(Updated_Seurat_Obj$New_Anno))
  subset_our_obj2 <- subset(Updated_Seurat_Obj,
                            cells = rownames(Updated_Seurat_Obj@meta.data)[which(Updated_Seurat_Obj$New_Anno %in% c("E18_HSPC", "E18_Stroma"))])
  
  ### run PCA and UMAP
  subset_our_obj2 <- FindVariableFeatures(subset_our_obj2)
  subset_our_obj2 <- ScaleData(subset_our_obj2)
  subset_our_obj2 <- RunPCA(subset_our_obj2, npcs = 15)
  ElbowPlot(subset_our_obj2, ndims = 15, reduction = "pca")
  subset_our_obj2 <- RunUMAP(subset_our_obj2, dims = 1:15)
  
  ### run RNAMagnet with anchors
  subset_our_obj2 <- SetIdent(object = subset_our_obj2,
                              cells = rownames(subset_our_obj2@meta.data),
                              value = subset_our_obj2$New_Anno)
  subset_our_obj2$New_Anno <- factor(subset_our_obj2$New_Anno,
                                     levels = as.character(unique(subset_our_obj2$New_Anno)))
  our_result_E18 <- RNAMagnetAnchors(seurat = subset_our_obj2,
                                     anchors = levels(subset_our_obj2$New_Anno),
                                     human = FALSE)
  
  ### write the result as an Excel file
  write.xlsx2(data.frame(Cell=rownames(our_result_E18), our_result_E18,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "RNAMagnet_Anchors_Ours_E18.xlsx"),
              sheetName = paste0("RNAMagnet_Result"),
              row.names = FALSE)
  
  ### add RNAMagnet info to the seurat object
  subset_our_obj2$direction <- as.character(our_result_E18[rownames(subset_our_obj2@meta.data),"direction"])
  subset_our_obj2$adhesiveness <- as.numeric(our_result_E18[rownames(subset_our_obj2@meta.data),"adhesiveness"])
  subset_our_obj2$specificity <- as.numeric(sapply(rownames(subset_our_obj2@meta.data), function(x) {
    our_result_E18[x, our_result_E18[x,"direction"]]
  }))
  
  ### make a data frame for ggplot
  dim_map <- Embeddings(subset_our_obj2, reduction = "umap")[rownames(subset_our_obj2@meta.data),]
  plot_df <- data.frame(X=dim_map[rownames(subset_our_obj2@meta.data),1],
                        Y=dim_map[rownames(subset_our_obj2@meta.data),2],
                        cluster_color = subset_our_obj2$New_Anno,
                        direction = subset_our_obj2@meta.data$direction,
                        adhesiveness = subset_our_obj2@meta.data$adhesiveness,
                        specificity = subset_our_obj2@meta.data$specificity,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### draw a scatter plot with the adhesiveness info
  p <- list()
  p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
    xlab("UMAP1") + ylab("UMAP2") +
    labs(col="Cluster") +
    ggtitle("UMAP with Cell Type") +
    scale_color_brewer(palette="Dark2") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'))
  ### the id column in the plot_df should be a factor
  # p[[1]] <- LabelClusters(plot = p[[1]], id = "cluster_color", col = "black")
  
  p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="direction", alpha="specificity"), size=2) +
    xlab("UMAP1") + ylab("UMAP2") +
    labs(col="Direction", alpha="Specificity") +
    ggtitle("UMAP with Direction & Specificity") +
    scale_color_brewer(palette="Set1") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'))
  ### the id column in the plot_df should be a factor
  # p[[2]] <- LabelClusters(plot = p[[3]], id = "cluster_color", col = "black")
  
  ### save the plots
  g <- arrangeGrob(grobs = p,
                   nrow = 1,
                   ncol = 2,
                   top = "")
  ggsave(file = paste0(outputDir, "RNAMagnet_Anchors_Ours_E18.png"), g, width = 30, height = 10, dpi = 350)
  
  
  ### pull out P0 HSC/MPP and Stroma only
  print(unique(Updated_Seurat_Obj$New_Anno))
  subset_our_obj3 <- subset(Updated_Seurat_Obj,
                            cells = rownames(Updated_Seurat_Obj@meta.data)[which(Updated_Seurat_Obj$New_Anno %in% c("P0_HSPC", "P0_Stroma"))])
  
  ### run PCA and UMAP
  subset_our_obj3 <- FindVariableFeatures(subset_our_obj3)
  subset_our_obj3 <- ScaleData(subset_our_obj3)
  subset_our_obj3 <- RunPCA(subset_our_obj3, npcs = 15)
  ElbowPlot(subset_our_obj3, ndims = 15, reduction = "pca")
  subset_our_obj3 <- RunUMAP(subset_our_obj3, dims = 1:15)
  
  ### run RNAMagnet with anchors
  subset_our_obj3 <- SetIdent(object = subset_our_obj3,
                              cells = rownames(subset_our_obj3@meta.data),
                              value = subset_our_obj3$New_Anno)
  subset_our_obj3$New_Anno <- factor(subset_our_obj3$New_Anno,
                                     levels = as.character(unique(subset_our_obj3$New_Anno)))
  our_result_P0 <- RNAMagnetAnchors(seurat = subset_our_obj3,
                                    anchors = levels(subset_our_obj3$New_Anno),
                                    human = FALSE)
  
  ### write the result as an Excel file
  write.xlsx2(data.frame(Cell=rownames(our_result_P0), our_result_P0,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "RNAMagnet_Anchors_Ours_P0.xlsx"),
              sheetName = paste0("RNAMagnet_Result"),
              row.names = FALSE)
  
  ### add RNAMagnet info to the seurat object
  subset_our_obj3$direction <- as.character(our_result_P0[rownames(subset_our_obj3@meta.data),"direction"])
  subset_our_obj3$adhesiveness <- as.numeric(our_result_P0[rownames(subset_our_obj3@meta.data),"adhesiveness"])
  subset_our_obj3$specificity <- as.numeric(sapply(rownames(subset_our_obj3@meta.data), function(x) {
    our_result_P0[x, our_result_P0[x,"direction"]]
  }))
  
  ### make a data frame for ggplot
  dim_map <- Embeddings(subset_our_obj3, reduction = "umap")[rownames(subset_our_obj3@meta.data),]
  plot_df <- data.frame(X=dim_map[rownames(subset_our_obj3@meta.data),1],
                        Y=dim_map[rownames(subset_our_obj3@meta.data),2],
                        cluster_color = subset_our_obj3$New_Anno,
                        direction = subset_our_obj3@meta.data$direction,
                        adhesiveness = subset_our_obj3@meta.data$adhesiveness,
                        specificity = subset_our_obj3@meta.data$specificity,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### draw a scatter plot with the adhesiveness info
  p <- list()
  p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
    xlab("UMAP1") + ylab("UMAP2") +
    labs(col="Cluster") +
    ggtitle("UMAP with Cell Type") +
    scale_color_brewer(palette="Dark2") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'))
  ### the id column in the plot_df should be a factor
  # p[[1]] <- LabelClusters(plot = p[[1]], id = "cluster_color", col = "black")
  
  p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="direction", alpha="specificity"), size=2) +
    xlab("UMAP1") + ylab("UMAP2") +
    labs(col="Direction", alpha="Specificity") +
    ggtitle("UMAP with Direction & Specificity") +
    scale_color_brewer(palette="Set1") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'))
  ### the id column in the plot_df should be a factor
  # p[[2]] <- LabelClusters(plot = p[[3]], id = "cluster_color", col = "black")
  
  ### save the plots
  g <- arrangeGrob(grobs = p,
                   nrow = 1,
                   ncol = 2,
                   top = "")
  ggsave(file = paste0(outputDir, "RNAMagnet_Anchors_Ours_P0.png"), g, width = 30, height = 10, dpi = 350)
  
  
  ### RNAMagnet signaling with our E16
  result_E16 <- RNAMagnetSignaling(seurat = subset_our_obj1,
                                   human = FALSE)
  
  ### get all the interaction list
  interaction_list_E16 <- NULL
  for(clust1 in levels(subset_our_obj1$New_Anno)) {
    for(clust2 in levels(subset_our_obj1$New_Anno)) {
      il <- getRNAMagnetGenes(result_E16, clust1, clust2, thresh = 0)
      if(nrow(il) > 0) {
        il$ligand_cluster <- clust1
        il$receptor_cluster <- clust2
        if(is.null(interaction_list_E16)) {
          interaction_list_E16 <- il
        } else {
          interaction_list_E16 <- rbind(interaction_list_E16, il)
        }
      }
    }
  }
  interaction_list_E16 <- cbind(interaction_list_E16[3:4], interaction_list_E16[1:2])
  colnames(interaction_list_E16) <- c("Ligand_Cluster", "Receptor_Cluster", "Interaction_Score", "Interaction_Pair")
  
  ### write out the interaction list
  write.xlsx2(interaction_list_E16,
              file = paste0(outputDir, "RNAMagnet_Signaling_Ours_E16.xlsx"),
              sheetName = paste0("Signaling_RNAMagnet_Interaction_List"),
              row.names = FALSE)
  
  
  ### RNAMagnet signaling with our E18
  result_E18 <- RNAMagnetSignaling(seurat = subset_our_obj2,
                                   human = FALSE)
  
  ### get all the interaction list
  interaction_list_E18 <- NULL
  for(clust1 in levels(subset_our_obj2$New_Anno)) {
    for(clust2 in levels(subset_our_obj2$New_Anno)) {
      il <- getRNAMagnetGenes(result_E18, clust1, clust2, thresh = 0)
      if(nrow(il) > 0) {
        il$ligand_cluster <- clust1
        il$receptor_cluster <- clust2
        if(is.null(interaction_list_E18)) {
          interaction_list_E18 <- il
        } else {
          interaction_list_E18 <- rbind(interaction_list_E18, il)
        }
      }
    }
  }
  interaction_list_E18 <- cbind(interaction_list_E18[3:4], interaction_list_E18[1:2])
  colnames(interaction_list_E18) <- c("Ligand_Cluster", "Receptor_Cluster", "Interaction_Score", "Interaction_Pair")
  
  ### write out the interaction list
  write.xlsx2(interaction_list_E18,
              file = paste0(outputDir, "RNAMagnet_Signaling_Ours_E18.xlsx"),
              sheetName = paste0("Signaling_RNAMagnet_Interaction_List"),
              row.names = FALSE)
  
  
  ### RNAMagnet signaling with our P0
  result_P0 <- RNAMagnetSignaling(seurat = subset_our_obj3,
                                  human = FALSE)
  
  ### get all the interaction list
  interaction_list_P0 <- NULL
  for(clust1 in levels(subset_our_obj3$New_Anno)) {
    for(clust2 in levels(subset_our_obj3$New_Anno)) {
      il <- getRNAMagnetGenes(result_P0, clust1, clust2, thresh = 0)
      if(nrow(il) > 0) {
        il$ligand_cluster <- clust1
        il$receptor_cluster <- clust2
        if(is.null(interaction_list_P0)) {
          interaction_list_P0 <- il
        } else {
          interaction_list_P0 <- rbind(interaction_list_P0, il)
        }
      }
    }
  }
  interaction_list_P0 <- cbind(interaction_list_P0[3:4], interaction_list_P0[1:2])
  colnames(interaction_list_P0) <- c("Ligand_Cluster", "Receptor_Cluster", "Interaction_Score", "Interaction_Pair")
  
  ### write out the interaction list
  write.xlsx2(interaction_list_P0,
              file = paste0(outputDir, "RNAMagnet_Signaling_Ours_P0.xlsx"),
              sheetName = paste0("Signaling_RNAMagnet_Interaction_List"),
              row.names = FALSE)
  
  
  ### correlation plot between the interaction scores
  ### ligand stroma, receptor MPP
  interaction_list2 <- interaction_list[intersect(which(interaction_list$Ligand_Cluster == "stroma"),
                                                  which(interaction_list$Receptor_Cluster == "HSC/MPP and pro")),]
  interaction_list2_E16 <- interaction_list_E16[intersect(which(interaction_list_E16$Ligand_Cluster == "E16_Stroma"),
                                                          which(interaction_list_E16$Receptor_Cluster == "E16_HSPC")),]
  interaction_list2_E18 <- interaction_list_E18[intersect(which(interaction_list_E18$Ligand_Cluster == "E18_Stroma"),
                                                          which(interaction_list_E18$Receptor_Cluster == "E18_HSPC")),]
  interaction_list2_P0 <- interaction_list_P0[intersect(which(interaction_list_P0$Ligand_Cluster == "P0_Stroma"),
                                                        which(interaction_list_P0$Receptor_Cluster == "P0_HSPC")),]
  
  temp1 <- merge(x = interaction_list2, y = interaction_list2_E16,
                 by.x = "Interaction_Pair", by.y = "Interaction_Pair",
                 all.x = FALSE, all.y = FALSE)
  temp1 <- temp1[,c("Interaction_Pair", "Interaction_Score.x", "Interaction_Score.y")]
  colnames(temp1)[2:3] <- c("Jardine", "Our_E16")
  temp2 <- merge(x = interaction_list2_E18, y = interaction_list2_P0,
                 by.x = "Interaction_Pair", by.y = "Interaction_Pair",
                 all.x = FALSE, all.y = FALSE)
  temp2 <- temp2[,c("Interaction_Pair", "Interaction_Score.x", "Interaction_Score.y")]
  colnames(temp2)[2:3] <- c("Our_E18", "Our_P0")
  
  plot_df <- merge(x = temp1, y = temp2,
                   by.x = "Interaction_Pair", by.y = "Interaction_Pair",
                   all.x = FALSE, all.y = FALSE)
  plot_df2 <- data.frame(Interaction_Pair=rep(plot_df$Interaction_Pair, 4),
                         Interaction_Score=c(plot_df$Jardine,
                                             plot_df$Our_E16,
                                             plot_df$Our_E18,
                                             plot_df$Our_P0),
                         Group=c(rep("Jardine", nrow(plot_df)),
                                 rep("E16", nrow(plot_df)),
                                 rep("E18", nrow(plot_df)),
                                 rep("P0", nrow(plot_df))),
                         stringsAsFactors = FALSE, check.names = FALSE)
  
  ### first, line plot of all the 4
  plot_df2$Group <- factor(plot_df2$Group, levels = c("E16", "E18", "P0", "Jardine"))
  p <- ggplot(data = plot_df2, aes_string(x="Interaction_Pair", y="Interaction_Score", group="Group")) +
    geom_point(aes_string(color="Group"), size = 1) +
    geom_line(aes_string(color="Group"), size = 2) +
    xlab("") +
    ylab("Interaction Score") +
    scale_color_brewer(palette="Dark2") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          axis.ticks = element_blank(),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'))
  ggsave(file = paste0(outputDir, "Interaction_Score_Comparison.png"), plot = p,
         width = 15, height = 10, dpi = 350)
  
  ### now pairwise correlation plots
  p <- list()
  ### E16
  s_cor <- round(cor(plot_df$Our_E16,
                     plot_df$Jardine, method = "spearman", use = "complete.obs"), 2)
  pv <- round(cor.test(plot_df$Our_E16,
                       plot_df$Jardine, method = "spearman", use = "complete.obs")$p.value, 2)
  p[[1]] <- ggplot(data = plot_df, aes(x=Our_E16, y=Jardine)) +
    geom_point(col = "aquamarine4", size = 5) +
    labs(title = paste0("Spearman Correlation:", s_cor),
         subtitle = paste0("P-value:", pv)) +
    xlab("Our E16") +
    ylab("Jardine") +
    geom_smooth(method = lm, color="coral3", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 40))
  ### E18
  s_cor <- round(cor(plot_df$Our_E18,
                     plot_df$Jardine, method = "spearman", use = "complete.obs"), 2)
  pv <- round(cor.test(plot_df$Our_E18,
                       plot_df$Jardine, method = "spearman", use = "complete.obs")$p.value, 2)
  p[[2]] <- ggplot(data = plot_df, aes(x=Our_E18, y=Jardine)) +
    geom_point(col = "aquamarine4", size = 5) +
    labs(title = paste0("Spearman Correlation:", s_cor),
         subtitle = paste0("P-value:", pv)) +
    xlab("Our E18") +
    ylab("Jardine") +
    geom_smooth(method = lm, color="coral3", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 40))
  ### P0
  s_cor <- round(cor(plot_df$Our_P0,
                     plot_df$Jardine, method = "spearman", use = "complete.obs"), 2)
  pv <- round(cor.test(plot_df$Our_P0,
                       plot_df$Jardine, method = "spearman", use = "complete.obs")$p.value, 2)
  p[[3]] <- ggplot(data = plot_df, aes(x=Our_P0, y=Jardine)) +
    geom_point(col = "aquamarine4", size = 5) +
    labs(title = paste0("Spearman Correlation:", s_cor),
         subtitle = paste0("P-value:", pv)) +
    xlab("Our P0") +
    ylab("Jardine") +
    geom_smooth(method = lm, color="coral3", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 40))
  
  ### save the plots
  g <- arrangeGrob(grobs = p,
                   nrow = 2,
                   ncol = 2,
                   top = "")
  ggsave(file = paste0(outputDir, "RNAMagnet_Signaling_Correlations.png"), g, width = 20, height = 15, dpi = 350)
  
  
  ### Pull out what they have annotated as HSC/MPP, recluster, and then perform GO analysis on the clusters
  ### to see if there is any overlap in the GO terms between their HSC/MPPs and our HSCs.
  subset_jardine_obj2 <- subset(jardine_seurat2,
                                cells = rownames(jardine_seurat2@meta.data)[which(jardine_seurat2$broad_fig1_cell.labels %in% c("HSC/MPP and pro"))])
  subset_jardine_obj2 <- FindVariableFeatures(subset_jardine_obj2)
  subset_jardine_obj2 <- ScaleData(subset_jardine_obj2)
  subset_jardine_obj2 <- RunPCA(subset_jardine_obj2, npcs = 15)
  ElbowPlot(subset_jardine_obj2, ndims = 15, reduction = "pca")
  subset_jardine_obj2 <- RunUMAP(subset_jardine_obj2, dims = 1:15)
  
  subset_jardine_obj2 <- FindNeighbors(subset_jardine_obj2, dims = 1:15)
  subset_jardine_obj2 <- FindClusters(subset_jardine_obj2, resolution = 0.2)
  
  p <- DimPlot(object = subset_jardine_obj2, reduction = "umap", raster = TRUE,
          group.by = "seurat_clusters",
          pt.size = 1) +
    ggtitle("") +
    labs(color = "seurat_clusters") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 48) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 36, color = "black", face = "bold"),
          axis.title = element_text(size = 36, color = "black", face = "bold"),
          axis.text = element_text(size = 36, color = "black", face = "bold"),
          legend.title = element_text(size = 24, color = "black", face = "bold"),
          legend.text = element_text(size = 24, color = "black", face = "bold"),
          axis.ticks = element_blank())
  ggsave(paste0(outputDir, "UMAP_Jardine_HSPCs.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  ### find all markers
  subset_jardine_obj2 <- SetIdent(object = subset_jardine_obj2,
                                  cells = rownames(subset_jardine_obj2@meta.data),
                                  value = subset_jardine_obj2$seurat_clusters)
  jardine_de_result <- FindAllMarkers(subset_jardine_obj2,
                                      min.pct = 0.2,
                                      logfc.threshold = 0.2,
                                      test.use = "wilcox")
  
  ### ours from E16, E18, and P0
  subset_our_obj4 <- subset(Updated_Seurat_Obj,
                            cells = rownames(Updated_Seurat_Obj@meta.data)[which(Updated_Seurat_Obj$New_HSPC %in% c("HSPC"))])
  
  subset_our_obj4 <- FindVariableFeatures(subset_our_obj4)
  subset_our_obj4 <- ScaleData(subset_our_obj4)
  subset_our_obj4 <- RunPCA(subset_our_obj4, npcs = 15)
  ElbowPlot(subset_our_obj4, ndims = 15, reduction = "pca")
  subset_our_obj4 <- RunUMAP(subset_our_obj4, dims = 1:15)
  
  subset_our_obj4 <- FindNeighbors(subset_our_obj4, dims = 1:15)
  subset_our_obj4 <- FindClusters(subset_our_obj4, resolution = 0.2)
  
  p <- DimPlot(object = subset_our_obj4, reduction = "umap", raster = TRUE,
               group.by = "seurat_clusters",
               pt.size = 1) +
    ggtitle("") +
    labs(color = "seurat_clusters") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 48) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 36, color = "black", face = "bold"),
          axis.title = element_text(size = 36, color = "black", face = "bold"),
          axis.text = element_text(size = 36, color = "black", face = "bold"),
          legend.title = element_text(size = 24, color = "black", face = "bold"),
          legend.text = element_text(size = 24, color = "black", face = "bold"),
          axis.ticks = element_blank())
  ggsave(paste0(outputDir, "UMAP_Ours_HSPCs.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  ### find all markers
  subset_our_obj4 <- SetIdent(object = subset_our_obj4,
                              cells = rownames(subset_our_obj4@meta.data),
                              value = subset_our_obj4$seurat_clusters)
  our_de_result <- FindAllMarkers(subset_our_obj4,
                                  min.pct = 0.2,
                                  logfc.threshold = 0.2,
                                  test.use = "wilcox")
  
  ### now compare DE genes between Jardine and Ours
  similarity_mat <- matrix(0, nrow = length(levels(jardine_de_result$cluster)), ncol = length(levels(our_de_result$cluster)))
  rownames(similarity_mat) <- paste0("Jardine_", levels(jardine_de_result$cluster))
  colnames(similarity_mat) <- paste0("Ours_", levels(our_de_result$cluster))
  
  ### calculate common DE genes - Jaccard Similarity (adjuated pv < 0.05)
  for(clstr1 in levels(jardine_de_result$cluster)) {
    for(clstr2 in levels(our_de_result$cluster)) {
      jardine_target_genes <- jardine_de_result$gene[intersect(which(jardine_de_result$p_val_adj < 0.05),
                                                               which(jardine_de_result$cluster == clstr1))]
      jardine_target_genes <- convertHumanGeneList(jardine_target_genes)
      our_target_genes <- our_de_result$gene[intersect(which(our_de_result$p_val_adj < 0.05),
                                                       which(our_de_result$cluster == clstr2))]
      
      similarity_mat[paste0("Jardine_", clstr1),paste0("Ours_", clstr2)] <- length(intersect(jardine_target_genes,
                                                                                             our_target_genes)) / length(union(jardine_target_genes,
                                                                                                                               our_target_genes))
    }
  }
  
  ### write out the interaction list
  write.xlsx2(data.frame(Cluster=rownames(similarity_mat),
                         similarity_mat,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "Cluster_DE_Gene_Comparison_JS.xlsx"),
              sheetName = paste0("Cluster_DE_Gene_Comparison"),
              row.names = FALSE)
  
  
  # ******************************************************************************************
  # Pathway Analysis with clusterProfiler package
  # Input: geneList     = a vector of gene Entrez IDs for pathway analysis [numeric or character]
  #        org          = organism that will be used in the analysis ["human" or "mouse"]
  #                       should be either "human" or "mouse"
  #        database     = pathway analysis database (KEGG or GO) ["KEGG" or "GO"]
  #        title        = title of the pathway figure [character]
  #        pv_threshold = pathway analysis p-value threshold (not DE analysis threshold) [numeric]
  #        displayNum   = the number of pathways that will be displayed [numeric]
  #                       (If there are many significant pathways show the few top pathways)
  #        imgPrint     = print a plot of pathway analysis [TRUE/FALSE]
  #        dir          = file directory path of the output pathway figure [character]
  #
  # Output: Pathway analysis results in figure - using KEGG and GO pathways
  #         The x-axis represents the number of DE genes in the pathway
  #         The y-axis represents pathway names
  #         The color of a bar indicates adjusted p-value from the pathway analysis
  #         For Pathview Result, all colored genes are found DE genes in the pathway,
  #         and the color indicates log2(fold change) of the DE gene from DE analysis
  # ******************************************************************************************
  pathwayAnalysis_CP <- function(geneList,
                                 org,
                                 database,
                                 title="Pathway_Results",
                                 pv_threshold=0.05,
                                 displayNum=Inf,
                                 imgPrint=TRUE,
                                 dir="./") {
    
    ### load library
    if(!require(clusterProfiler, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("clusterProfiler")
      require(clusterProfiler, quietly = TRUE)
    }
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    
    
    ### collect gene list (Entrez IDs)
    geneList <- geneList[which(!is.na(geneList))]
    
    if(!is.null(geneList)) {
      ### make an empty list
      p <- list()
      
      if(database == "KEGG") {
        ### KEGG Pathway
        kegg_enrich <- enrichKEGG(gene = geneList, organism = org, pvalueCutoff = pv_threshold)
        
        if(is.null(kegg_enrich)) {
          writeLines("KEGG Result does not exist")
          return(NULL)
        } else {
          kegg_enrich@result <- kegg_enrich@result[which(kegg_enrich@result$p.adjust < pv_threshold),]
          
          if(imgPrint == TRUE) {
            if((displayNum == Inf) || (nrow(kegg_enrich@result) <= displayNum)) {
              result <- kegg_enrich@result
              description <- kegg_enrich@result$Description
            } else {
              result <- kegg_enrich@result[1:displayNum,]
              description <- kegg_enrich@result$Description[1:displayNum]
            }
            
            if(nrow(kegg_enrich) > 0) {
              p[[1]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
                theme_classic(base_size = 30) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("KEGG ", title)) +
                theme(axis.text = element_text(size = 30))
              
              ggsave(file = paste0(dir, "KEGG_", title, "_CB.png"), plot = p[[1]], width = 35, height = 10, dpi = 350)
            } else {
              writeLines("KEGG Result does not exist")
            }
          }
          
          return(kegg_enrich@result)
        }
      } else if(database == "GO") {
        ### GO Pathway
        if(org == "human") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Hs.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else if(org == "mouse") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Mm.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else {
          go_enrich <- NULL
          writeLines(paste("Unknown org variable:", org))
        }
        
        if(is.null(go_enrich)) {
          writeLines("GO Result does not exist")
          return(NULL)
        } else {
          go_enrich@result <- go_enrich@result[which(go_enrich@result$p.adjust < pv_threshold),]
          
          if(imgPrint == TRUE) {
            if((displayNum == Inf) || (nrow(go_enrich@result) <= displayNum)) {
              result <- go_enrich@result
              description <- go_enrich@result$Description
            } else {
              result <- go_enrich@result[1:displayNum,]
              description <- go_enrich@result$Description[1:displayNum]
            }
            
            if(nrow(go_enrich) > 0) {
              p[[2]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
                theme_classic(base_size = 30) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("GO ", title)) +
                theme(axis.text = element_text(size = 30))
              
              ggsave(file = paste0(dir, "GO_", title, "_CB.png"), plot = p[[2]], width = 35, height = 10, dpi = 350)
            } else {
              writeLines("GO Result does not exist")
            }
          }
          
          return(go_enrich@result)
        }
      } else {
        stop("database prameter should be \"GO\" or \"KEGG\"")
      }
    } else {
      writeLines("geneList = NULL")
    }
  }
  
  
  ### now let's see the common GO terms
  Jardine_Up_Go <- vector("list", length = length(levels(jardine_de_result$cluster)))
  names(Jardine_Up_Go) <- levels(jardine_de_result$cluster)
  Jardine_Down_Go <- vector("list", length = length(levels(jardine_de_result$cluster)))
  names(Jardine_Down_Go) <- levels(jardine_de_result$cluster)
  Our_Up_Go <- vector("list", length = length(levels(our_de_result$cluster)))
  names(Our_Up_Go) <- levels(our_de_result$cluster)
  Our_Down_Go <- vector("list", length = length(levels(our_de_result$cluster)))
  names(Our_Down_Go) <- levels(our_de_result$cluster)
  
  ### Jardine
  for(clstr in levels(jardine_de_result$cluster)) {
    ### Up-Regulated
    jardine_target_up_genes <- jardine_de_result$gene[intersect(intersect(which(jardine_de_result$p_val_adj < 0.05),
                                                                          which(jardine_de_result$avg_log2FC > 0)),
                                                                which(jardine_de_result$cluster == clstr))]
    
    ### get entrez ids for the genes
    de_entrez_ids <- mapIds(org.Hs.eg.db,
                            jardine_target_up_genes,
                            "ENTREZID", "SYMBOL")
    de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
    
    Jardine_Up_Go[[clstr]] <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                                 org = "human", database = "GO",
                                                 imgPrint = FALSE)
    
    ### Down-Regulated
    jardine_target_down_genes <- jardine_de_result$gene[intersect(intersect(which(jardine_de_result$p_val_adj < 0.05),
                                                                            which(jardine_de_result$avg_log2FC < 0)),
                                                                  which(jardine_de_result$cluster == clstr))]
    
    ### get entrez ids for the genes
    de_entrez_ids <- mapIds(org.Hs.eg.db,
                            jardine_target_down_genes,
                            "ENTREZID", "SYMBOL")
    de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
    
    Jardine_Down_Go[[clstr]] <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                                   org = "human", database = "GO",
                                                   imgPrint = FALSE)
  }
  
  ### Ours
  for(clstr in levels(our_de_result$cluster)) {
    ### Up-Regulated
    our_target_up_genes <- our_de_result$gene[intersect(intersect(which(our_de_result$p_val_adj < 0.05),
                                                                  which(our_de_result$avg_log2FC > 0)),
                                                        which(our_de_result$cluster == clstr))]
    
    ### get entrez ids for the genes
    de_entrez_ids <- mapIds(org.Mm.eg.db,
                            our_target_up_genes,
                            "ENTREZID", "SYMBOL")
    de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
    
    Our_Up_Go[[clstr]] <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                             org = "mouse", database = "GO",
                                             imgPrint = FALSE)
    
    ### Down-Regulated
    our_target_down_genes <- our_de_result$gene[intersect(intersect(which(our_de_result$p_val_adj < 0.05),
                                                                    which(our_de_result$avg_log2FC < 0)),
                                                          which(our_de_result$cluster == clstr))]
    
    ### get entrez ids for the genes
    de_entrez_ids <- mapIds(org.Mm.eg.db,
                            our_target_down_genes,
                            "ENTREZID", "SYMBOL")
    de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
    
    Our_Down_Go[[clstr]] <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                               org = "mouse", database = "GO",
                                               imgPrint = FALSE)
  }
  
  ### Pairwise GO term comparison
  
  
  
  
  
}

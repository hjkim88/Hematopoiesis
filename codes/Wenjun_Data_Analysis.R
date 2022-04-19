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

wenjun_analysis <- function(Seurat_Obj_Path="./data/Wenjun_Seurat_Obj2.RDS",
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
  ggsave(paste0(outputDir, "/Wenjun_UMAP_singleR_celltypes_Blueprint2(2).png"), plot = p, width = 20, height = 10, dpi = 350)
  
  ### CXCL12, FOXC1, EBF3 UMAP
  p <- FeaturePlot(WJ_Seurat_Obj, features = c("Cxcl12", "Foxc1", "Ebf3"), cols = c("lightgray", "red"))
  ggsave(paste0(outputDir, "/Wenjun_FeaturePlot_CAR_cells(2).png"), plot = p, width = 18, height = 10, dpi = 400)
  
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
  
  ### annotation by clusters
  WJ_Seurat_Obj$WJ_Annotation <- NA
  WJ_Seurat_Obj$WJ_Annotation[which(WJ_Seurat_Obj$seurat_clusters == "0")] <- "Chondrocytes"
  WJ_Seurat_Obj$WJ_Annotation[which(WJ_Seurat_Obj$seurat_clusters == "1")] <- "Monocyte progenitors"
  WJ_Seurat_Obj$WJ_Annotation[which(WJ_Seurat_Obj$seurat_clusters == "2")] <- "Fibroblasts"
  WJ_Seurat_Obj$WJ_Annotation[which(WJ_Seurat_Obj$seurat_clusters == "3")] <- "Neutrophils"
  WJ_Seurat_Obj$WJ_Annotation[which(WJ_Seurat_Obj$seurat_clusters == "4")] <- "Erythroblasts"
  WJ_Seurat_Obj$WJ_Annotation[which(WJ_Seurat_Obj$seurat_clusters == "5")] <- "B cells"
  WJ_Seurat_Obj$WJ_Annotation[which(WJ_Seurat_Obj$seurat_clusters == "6")] <- "B cells"
  WJ_Seurat_Obj$WJ_Annotation[which(WJ_Seurat_Obj$seurat_clusters == "7")] <- "Endothelial cells"
  WJ_Seurat_Obj$WJ_Annotation[which(WJ_Seurat_Obj$seurat_clusters == "8")] <- "Schwann cells"
  WJ_Seurat_Obj$WJ_Annotation[which(WJ_Seurat_Obj$seurat_clusters == "9")] <- "Smooth muscle cells"
  WJ_Seurat_Obj$WJ_Annotation[which(WJ_Seurat_Obj$seurat_clusters == "10")] <- "Myofibroblasts"
  WJ_Seurat_Obj$WJ_Annotation[which(WJ_Seurat_Obj$seurat_clusters == "11")] <- "CAR"
  WJ_Seurat_Obj$WJ_Annotation[which(WJ_Seurat_Obj$seurat_clusters == "12")] <- "Monocyte"
  
  cluster_anno <- c("Chondrocytes", "Monocyte progenitors", "Fibroblasts", "Neutrophils",
                    "Erythroblasts", "B cells", "B cells", "Endothelial cells", "Schwann cells",
                    "Smooth muscle cells", "Myofibroblasts", "CAR", "Monocyte")
  names(cluster_anno) <- 0:12
  
  p <- DimPlot(object = WJ_Seurat_Obj, reduction = "umap", raster = TRUE,
               group.by = "WJ_Annotation",
               pt.size = 1) +
    ggtitle(paste0("WJ Annotation")) +
    labs(color = "Annotation") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 48) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 36, color = "black", face = "bold"),
          axis.title = element_text(size = 36, color = "black", face = "bold"),
          axis.text = element_text(size = 36, color = "black", face = "bold"),
          legend.title = element_text(size = 24, color = "black", face = "bold"),
          legend.text = element_text(size = 24, color = "black", face = "bold"),
          axis.ticks = element_blank())
  ggsave(paste0(outputDir, "/Wenjun_UMAP_WJ_Annotation.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  
  ### Proportional Bar plot
  ### In each cluster, the proportion of different time point
  
  ### Create a data for the plot
  target_clusters <- levels(WJ_Seurat_Obj$seurat_clusters)
  target_times <- unique(WJ_Seurat_Obj$Development)
  target_times <- c("P6", "P14", "P21", "P28")
  plot_df <- data.frame(Cluster=as.vector(sapply(target_clusters, function(x) rep(x, length(target_times)))),
                        Time=rep(target_times, length(target_clusters)),
                        Num=0,
                        Pcnt=0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### fill out the data frame
  for(i in 1:nrow(plot_df)) {
    plot_df$Num[i] <- length(intersect(which(WJ_Seurat_Obj$seurat_clusters == plot_df$Cluster[i]),
                                       which(WJ_Seurat_Obj$Development == plot_df$Time[i])))
  }
  
  ### calculate percentages
  cluster_sum <- rep(0, length(target_clusters))
  names(cluster_sum) <- target_clusters
  for(i in 1:length(target_clusters)) {
    cluster_sum[i] <- sum(plot_df[which(plot_df$Cluster == target_clusters[i]),"Num"])
    plot_df$Pcnt[which(plot_df$Cluster == target_clusters[i])] <- round(plot_df$Num[which(plot_df$Cluster == target_clusters[i])] * 100 / cluster_sum[i], 1)
  }
  
  ### remove NaN rows
  plot_df <- plot_df[which(!is.nan(plot_df$Pcnt)),]
  
  ### preserve the pcnt
  plot_df$Pcnt2 <- plot_df$Pcnt
  
  ### pcnt < 5 -> ""
  plot_df$Pcnt[which(as.numeric(plot_df$Pcnt) < 5)] <- ""
  plot_df$Pcnt <- as.character(plot_df$Pcnt)
  
  ### add states2 column
  plot_df$Time2 <- ""
  plot_df$Time2[which(as.numeric(plot_df$Num) != 0)] <- as.character(plot_df$Time[which(as.numeric(plot_df$Num) != 0)])
  plot_df$Time2[which(plot_df$Pcnt == "")] <- ""
  
  ### annotate "%"
  plot_df$Pcnt[which(as.numeric(plot_df$Pcnt) != 0)] <- paste0(plot_df$Pcnt[which(as.numeric(plot_df$Pcnt) != 0)], "%")
  plot_df$Pcnt[which(plot_df$Pcnt == 0)] <- ""
  
  ### factorize the time point & state
  plot_df$Cluster <- factor(plot_df$Cluster, levels = target_clusters)
  plot_df$Time <- factor(plot_df$Time, levels = target_times)
  
  ### draw a proportional bar plot
  p <- ggplot(data=plot_df, aes_string(x="Cluster", y="Pcnt2", fill="Time", label="Pcnt")) +
    geom_bar(position = "stack", stat = "identity") +
    # ggtitle("Proportion of Cells") +
    xlab("Clusters") + ylab("Percentage") +
    # geom_text(size = 3, position = position_stack(vjust = 0.5), vjust = 3, color = "blue") +
    geom_text(aes_string(x="Cluster", y="Pcnt2", label = "Time2"),
              position = position_stack(vjust = 0.5),
              size = 12, color = "white") +
    coord_flip() +
    scale_x_discrete(labels = paste0(names(cluster_anno), " (", cluster_anno, ")")) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(size = 35, color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(0.7, 'cm'),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir, "Proportional_Barplot_Cluster_Time.png"), plot = p,
         width = 22, height = 10, dpi = 350)
  
  ### Cellchat
  
  ### create a cell chat object
  cellchat <- createCellChat(object = WJ_Seurat_Obj, group.by = "WJ_Annotation")
  
  ### set ligand-receptor interaction db
  ### CellChatDB <- CellChatDB.human
  CellChatDB <- CellChatDB.mouse
  
  ### use a subset of CellChatDB for cell-cell communication analysis
  ### CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  ### use all CellChatDB for cell-cell communication analysis
  CellChatDB.use <- CellChatDB
  
  ### set the used database in the object
  cellchat@DB <- CellChatDB.use
  
  ### subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  
  ### Preprocessing
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  ### project gene expression data onto PPI network (optional)
  cellchat <- projectData(cellchat, PPI.mouse)
  
  ### Compute the communication probability and infer cellular communication network
  cellchat <- computeCommunProb(cellchat)
  ### Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  ### Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  
  ### Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  
  ### We can also visualize the aggregated cell-cell communication network.
  ### For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
  groupSize <- as.numeric(table(cellchat@idents))
  png(filename = paste0(outputDir, "/", "Cellchat_Circle_Interaction.png"), width = 3000, height = 2500, res = 350)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  
  ### Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group.
  ### Here we also control the parameter edge.weight.max so that we can compare edge weights between different networks.
  mat <- cellchat@net$weight
  png(filename = paste0(outputDir, "/", "Cellchat_Circle_Interaction_Each.png"), width = 3200, height = 2500, res = 350)
  par(mfrow = c(3,4), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()
  
  ### One can visualize the inferred communication network of signaling pathways using netVisual_aggregate,
  ### and visualize the inferred communication networks of individual L-R pairs associated with that signaling pathway using netVisual_individual.
  print(cellchat@netP$pathways)
  pathways.show <- c("CXCL")
  png(filename = paste0(outputDir, "/", "Cellchat_Circle_CXCL.png"), width = 3000, height = 2500, res = 350)
  netVisual_aggregate(cellchat, signaling = pathways.show,  layout = "circle")
  dev.off()
  png(filename = paste0(outputDir, "/", "Cellchat_Heatmap_CXCL.png"), width = 3000, height = 2500, res = 350)
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  dev.off()
  
  ### Compute the contribution of each ligand-receptor pair to the overall signaling pathway and
  ### visualize cell-cell communication mediated by a single ligand-receptor pair
  png(filename = paste0(outputDir, "/", "Cellchat_LR_Contribution_CXCL.png"), width = 3000, height = 2500, res = 350)
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  dev.off()
  
  ### show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use')
  ### to other cell groups (defined by 'targets.use')
  png(filename = paste0(outputDir, "/", "Cellchat_Significant_LR.png"), width = 3000, height = 7000, res = 350)
  netVisual_bubble(cellchat, sources.use = 2, targets.use = c(3:6, 9, 11, 12), remove.isolate = TRUE, thresh = 0.01)
  dev.off()
  
  ### Gene expression of the signaling pathway genes
  png(filename = paste0(outputDir, "/", "Cellchat_Gene_Exp_CXCL.png"), width = 3000, height = 2500, res = 350)
  plotGeneExpression(cellchat, signaling = "CXCL")
  dev.off()
  
  ### RNA Magnet
  
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
  
  ### run RNAMagnet
  result <- RNAMagnetAnchors(seurat = WJ_Seurat_Obj,
                             anchors = unique(WJ_Seurat_Obj$WJ_Annotation),
                             human = FALSE)
  
  ### write the result as an Excel file
  write.xlsx2(data.frame(Cell=rownames(result), result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "RNAMagnet_Anchors_WJ.xlsx"),
              sheetName = paste0("RNAMagnet_Result"),
              row.names = FALSE)
  
  ### add RNAMagnet info to the seurat object
  WJ_Seurat_Obj$direction <- as.character(result[rownames(WJ_Seurat_Obj@meta.data),"direction"])
  WJ_Seurat_Obj$adhesiveness <- as.numeric(result[rownames(WJ_Seurat_Obj@meta.data),"adhesiveness"])
  WJ_Seurat_Obj$specificity <- as.numeric(sapply(rownames(WJ_Seurat_Obj@meta.data), function(x) {
    result[x, result[x,"direction"]]
  }))
  
  ### make a data frame for ggplot
  dim_map <- Embeddings(WJ_Seurat_Obj, reduction = "umap")[rownames(WJ_Seurat_Obj@meta.data),]
  plot_df <- data.frame(X=dim_map[rownames(WJ_Seurat_Obj@meta.data),1],
                        Y=dim_map[rownames(WJ_Seurat_Obj@meta.data),2],
                        cluster_color = WJ_Seurat_Obj$WJ_Annotation,
                        direction = WJ_Seurat_Obj$direction,
                        adhesiveness = WJ_Seurat_Obj$adhesiveness,
                        specificity = WJ_Seurat_Obj$specificity,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### RNAMagnet anchors UMAP plot
  p <- ggplot(plot_df, aes_string(x="X", y="Y")) +
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
  ggsave(file = paste0(outputDir, "RNAMagnet_Anchors_WJ.png"), p, width = 20, height = 10, dpi = 350)
  
  ### RNAMagnet signaling
  result2 <- RNAMagnetSignaling(seurat = WJ_Seurat_Obj,
                                human = FALSE)
  
  ### get all the interaction list
  interaction_list <- NULL
  for(clust1 in unique(WJ_Seurat_Obj$WJ_Annotation)) {
    for(clust2 in unique(WJ_Seurat_Obj$WJ_Annotation)) {
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
              file = paste0(outputDir, "RNAMagnet_Signaling_WJ.xlsx"),
              sheetName = paste0("Signaling_RNAMagnet_Interaction_List"),
              row.names = FALSE)
  
  
  ### add PO to the data
  previous_hspc_seurat <- readRDS(file = "./data/Combined_Seurat_Obj.RDS")
  
  ### check whether the orders are the same
  print(identical(rownames(previous_hspc_seurat@meta.data), colnames(previous_hspc_seurat@assays$RNA@counts)))
  print(identical(names(Idents(object = previous_hspc_seurat)), rownames(previous_hspc_seurat@meta.data)))
  
  ### create new column for the analysis
  previous_hspc_seurat@meta.data$Dev_Anno <- paste0(previous_hspc_seurat@meta.data$Development, "_",
                                                    previous_hspc_seurat@meta.data$Annotation)
  
  ### P0 object
  previous_hspc_seurat_po <- subset(previous_hspc_seurat,
                                    cells = rownames(previous_hspc_seurat@meta.data)[which(previous_hspc_seurat$Development == "P0")])
  
  ### combine
  WJ_Seurat_Obj <- merge(WJ_Seurat_Obj, previous_hspc_seurat_po)
  WJ_Seurat_Obj$Cell1[which(is.na(WJ_Seurat_Obj$Cell1))] <- WJ_Seurat_Obj$Cell_Type[which(is.na(WJ_Seurat_Obj$Cell1))]
  WJ_Seurat_Obj$Cell2[which(is.na(WJ_Seurat_Obj$Cell2))] <- paste0(WJ_Seurat_Obj$Development[which(is.na(WJ_Seurat_Obj$Cell2))],
                                                                         "_",
                                                                         WJ_Seurat_Obj$Cell1[which(is.na(WJ_Seurat_Obj$Cell2))])
  
  ### preprocess
  
  ### Basic function to convert human to mouse gene names
  convertHumanGeneList <- function(x){
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
    humanx <- unique(genesV2[, 2])
    
    return(humanx)
  }
  
  ### MT percentage (Human: ^MT-, Mouse: ^mt- or ^Mt-)
  WJ_Seurat_Obj[["percent.mt"]] <- PercentageFeatureSet(WJ_Seurat_Obj, pattern = "^mt-")
  
  ### Cell cycle score (will be used later for regression out)
  WJ_Seurat_Obj <- CellCycleScoring(object = WJ_Seurat_Obj,
                                    g2m.features = convertHumanGeneList(cc.genes$g2m.genes),
                                    s.features = convertHumanGeneList(cc.genes$s.genes))
  
  ### find variable genes
  WJ_Seurat_Obj <- FindVariableFeatures(WJ_Seurat_Obj,
                                        selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  WJ_Seurat_Obj <- ScaleData(WJ_Seurat_Obj,
                             vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
  
  ### PCA
  WJ_Seurat_Obj <- RunPCA(WJ_Seurat_Obj,
                          features = VariableFeatures(object = WJ_Seurat_Obj), npcs = 20)
  ElbowPlot(WJ_Seurat_Obj, ndims = 20, reduction = "pca")
  
  ### UMAP
  WJ_Seurat_Obj <- RunUMAP(WJ_Seurat_Obj, dims = 1:20)
  
  ### clustering
  WJ_Seurat_Obj <- FindNeighbors(WJ_Seurat_Obj, dims = 1:20)
  WJ_Seurat_Obj <- FindClusters(WJ_Seurat_Obj, resolution = 0.1)
  
  ### display cell types and time points together
  p <- DimPlot(object = WJ_Seurat_Obj, reduction = "umap", raster = TRUE,
               group.by = "Cell2",
               pt.size = 0.5) +
    ggtitle(paste0("")) +
    labs(color = "Cell Type") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 48) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 36, color = "black", face = "bold"),
          axis.title = element_text(size = 36, color = "black", face = "bold"),
          axis.text = element_text(size = 36, color = "black", face = "bold"),
          legend.title = element_text(size = 24, color = "black", face = "bold"),
          legend.text = element_text(size = 24, color = "black", face = "bold"),
          axis.ticks = element_blank())
  ggsave(paste0(outputDir, "/Wenjun_UMAP_Cell2_P0_Included.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  ### devide the object into different time points
  WJ_Seurat_Obj_TimeList <- vector("list", length(unique(WJ_Seurat_Obj$Development)))
  names(WJ_Seurat_Obj_TimeList) <- unique(WJ_Seurat_Obj$Development)
  
  for(time in unique(WJ_Seurat_Obj$Development)) {
    WJ_Seurat_Obj_TimeList[[time]] <- subset(WJ_Seurat_Obj,
                                             cells = rownames(WJ_Seurat_Obj@meta.data)[which(WJ_Seurat_Obj$Development == time)])
    
    ### cell chat
    ### create a cell chat object
    cellchat <- createCellChat(object = WJ_Seurat_Obj_TimeList[[time]], group.by = "WJ_Annotation")
    
    ### set ligand-receptor interaction db
    ### CellChatDB <- CellChatDB.human
    CellChatDB <- CellChatDB.mouse
    
    ### use a subset of CellChatDB for cell-cell communication analysis
    ### CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
    ### use all CellChatDB for cell-cell communication analysis
    CellChatDB.use <- CellChatDB
    
    ### set the used database in the object
    cellchat@DB <- CellChatDB.use
    
    ### subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    
    ### Preprocessing
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    ### project gene expression data onto PPI network (optional)
    cellchat <- projectData(cellchat, PPI.mouse)
    
    ### Compute the communication probability and infer cellular communication network
    cellchat <- computeCommunProb(cellchat)
    ### Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    ### Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)
    
    ### Calculate the aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    
    
  }
  
  
  ###
  ### 04/01/2022
  ### divide the seurat into Heme and Stroma
  Heme_WJ_Obj <- subset(WJ_Seurat_Obj,
                        cells = rownames(WJ_Seurat_Obj@meta.data)[which(WJ_Seurat_Obj$Cell1 == "Heme")])
  Stroma_WJ_Obj <- subset(WJ_Seurat_Obj,
                          cells = rownames(WJ_Seurat_Obj@meta.data)[which(WJ_Seurat_Obj$Cell1 == "Stroma")])
  
  ### run FindAllMarkers() in each cell type
  ### Heme
  DimPlot(object = Heme_WJ_Obj, reduction = "umap", raster = TRUE,
          group.by = "seurat_clusters",
          pt.size = 1, label = TRUE)
  sapply(levels(WJ_Seurat_Obj$seurat_clusters), function(x) {
    return(length(intersect(which(WJ_Seurat_Obj$seurat_clusters == x),
                            which(WJ_Seurat_Obj$Cell1 == "Heme"))))
  })
  Heme_WJ_Obj <- SetIdent(object = Heme_WJ_Obj,
                          cells = rownames(Heme_WJ_Obj@meta.data),
                          value = Heme_WJ_Obj$seurat_clusters)
  heme_de_result_all <- FindAllMarkers(Heme_WJ_Obj,
                                       min.cells.group = 10,
                                       min.pct = 0.2,
                                       logfc.threshold = 0.2,
                                       test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(heme_de_result_all),
                         heme_de_result_all,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/Wenjun_Cluster_AllMarkers_Heme.xlsx"),
              sheetName = "Wenjun_Cluster_AllMarkers_Heme", row.names = FALSE)
  
  # heme_de_result_all <- read.xlsx2(file = paste0(outputDir, "/Wenjun_Cluster_AllMarkers_Heme.xlsx"), sheetIndex = 1,
  #                                  stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
  
  
  ### Stroma
  DimPlot(object = Stroma_WJ_Obj, reduction = "umap", raster = TRUE,
          group.by = "seurat_clusters",
          pt.size = 1, label = TRUE)
  sapply(levels(WJ_Seurat_Obj$seurat_clusters), function(x) {
    return(length(intersect(which(WJ_Seurat_Obj$seurat_clusters == x),
                            which(WJ_Seurat_Obj$Cell1 == "Stroma"))))
  })
  Stroma_WJ_Obj <- SetIdent(object = Stroma_WJ_Obj,
                            cells = rownames(Stroma_WJ_Obj@meta.data),
                            value = Stroma_WJ_Obj$seurat_clusters)
  stroma_de_result_all <- FindAllMarkers(Stroma_WJ_Obj,
                                         min.cells.group = 20,
                                         min.pct = 0.2,
                                         logfc.threshold = 0.2,
                                         test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(stroma_de_result_all),
                         stroma_de_result_all,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/Wenjun_Cluster_AllMarkers_Stroma.xlsx"),
              sheetName = "Wenjun_Cluster_AllMarkers_Stroma", row.names = FALSE)
  
  # stroma_de_result_all <- read.xlsx2(file = paste0(outputDir, "/Wenjun_Cluster_AllMarkers_Stroma.xlsx"), sheetIndex = 1,
  #                                    stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
  
  
  ### put min.cells.group worked well?
  ### Heme should not contain cluster 6, 13, 20, 23
  ### Stroma should not contain cluster 3, 4, 5, 9, 16, 18, 22
  print(unique(heme_de_result_all$cluster))
  print(unique(stroma_de_result_all$cluster))
  
  
  ###
  ### computational cluster annotation
  ###
  
  ### load reference annotation file
  ref_anno <- read.xlsx2(file = paste0("./data/Chiara_Baccin_NCB_2019_marker_gene.xlsx"),
                         sheetIndex = 1,
                         stringsAsFactors = FALSE, check.names = FALSE)
  ref_anno_list <- apply(ref_anno, 2, function(x) {
    y <- x
    y <- y[which(y != "")]
    return(y)
  })
  
  ### jacard index
  cal_jaccard <- function(x, y) {
    return(length(intersect(x, y)) / length(union(x, y)))
  }
  
  ### see the shared markers
  ### calculate jaccard index between clusters and cell types
  heme_jaccard_mat <- matrix(0, nrow = length(unique(heme_de_result_all$cluster)), ncol = length(ref_anno_list))
  rownames(heme_jaccard_mat) <- unique(heme_de_result_all$cluster)
  colnames(heme_jaccard_mat) <- names(ref_anno_list)
  heme_pos_jaccard_mat <- matrix(0, nrow = length(unique(heme_de_result_all$cluster)), ncol = length(ref_anno_list))
  rownames(heme_pos_jaccard_mat) <- unique(heme_de_result_all$cluster)
  colnames(heme_pos_jaccard_mat) <- names(ref_anno_list)
  for(clstr in rownames(heme_jaccard_mat)) {
    for(cell_type in colnames(heme_jaccard_mat)) {
      clstr_markers <- heme_de_result_all$gene[intersect(which(heme_de_result_all$cluster == clstr),
                                                         which(heme_de_result_all$p_val_adj < 0.05))]
      clstr_pos_markers <- intersect(clstr_markers,
                                     heme_de_result_all$gene[which(heme_de_result_all$avg_log2FC > 0)])
      cell_type_markers <- ref_anno_list[[cell_type]]
      heme_jaccard_mat[clstr,cell_type] <- cal_jaccard(clstr_markers, cell_type_markers)
      heme_pos_jaccard_mat[clstr,cell_type] <- cal_jaccard(clstr_pos_markers, cell_type_markers)
    }
  }
  
  stroma_jaccard_mat <- matrix(0, nrow = length(unique(stroma_de_result_all$cluster)), ncol = length(ref_anno_list))
  rownames(stroma_jaccard_mat) <- unique(stroma_de_result_all$cluster)
  colnames(stroma_jaccard_mat) <- names(ref_anno_list)
  stroma_pos_jaccard_mat <- matrix(0, nrow = length(unique(stroma_de_result_all$cluster)), ncol = length(ref_anno_list))
  rownames(stroma_pos_jaccard_mat) <- unique(stroma_de_result_all$cluster)
  colnames(stroma_pos_jaccard_mat) <- names(ref_anno_list)
  for(clstr in rownames(stroma_jaccard_mat)) {
    for(cell_type in colnames(stroma_jaccard_mat)) {
      clstr_markers <- stroma_de_result_all$gene[intersect(which(stroma_de_result_all$cluster == clstr),
                                                           which(stroma_de_result_all$p_val_adj < 0.05))]
      clstr_pos_markers <- intersect(clstr_markers,
                                     stroma_de_result_all$gene[which(stroma_de_result_all$avg_log2FC > 0)])
      cell_type_markers <- ref_anno_list[[cell_type]]
      stroma_jaccard_mat[clstr,cell_type] <- cal_jaccard(clstr_markers, cell_type_markers)
      stroma_pos_jaccard_mat[clstr,cell_type] <- cal_jaccard(clstr_pos_markers, cell_type_markers)
    }
  }
  
  ### first, second, third guess
  heme_guessing <- sapply(rownames(heme_jaccard_mat), function(x) {
    temp_row <- heme_jaccard_mat[x,]
    temp_row <- temp_row[order(-temp_row)]
    return(names(temp_row)[1:3])
  })
  
  heme_pos_guessing <- sapply(rownames(heme_pos_jaccard_mat), function(x) {
    temp_row <- heme_pos_jaccard_mat[x,]
    temp_row <- temp_row[order(-temp_row)]
    return(names(temp_row)[1:3])
  })
  
  stroma_guessing <- sapply(rownames(stroma_jaccard_mat), function(x) {
    temp_row <- stroma_jaccard_mat[x,]
    temp_row <- temp_row[order(-temp_row)]
    return(names(temp_row)[1:3])
  })
  
  stroma_pos_guessing <- sapply(rownames(stroma_pos_jaccard_mat), function(x) {
    temp_row <- stroma_pos_jaccard_mat[x,]
    temp_row <- temp_row[order(-temp_row)]
    return(names(temp_row)[1:3])
  })
  
  ### check how the normal and the pos different from each other
  ### it turned out they are exactly the same
  print(identical(heme_guessing, heme_pos_guessing))
  print(identical(stroma_guessing, stroma_pos_guessing))
  
  ### computational annotation
  ### cluster 7 is shared but should be Heme cells
  ### Anyway both Heme & Stroma cluster 7 indicates the same cell type, so no problem
  WJ_Seurat_Obj$computational_annotation <- NA
  for(clstr in colnames(stroma_guessing)) {
    WJ_Seurat_Obj$computational_annotation[which(WJ_Seurat_Obj$seurat_clusters == clstr)] <- stroma_guessing[1,clstr]
  }
  for(clstr in colnames(heme_guessing)) {
    WJ_Seurat_Obj$computational_annotation[which(WJ_Seurat_Obj$seurat_clusters == clstr)] <- heme_guessing[1,clstr]
  }
  
  ### UMAP plot with new computational annotation
  p <- DimPlot(object = WJ_Seurat_Obj, reduction = "umap", raster = TRUE,
               group.by = "computational_annotation",
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
  ggsave(paste0(outputDir, "/Wenjun_UMAP_Computational_Annotation.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  ###
  ### ok. now, cellchat for each time point
  ###
  
  cellchat_list <- vector("list", length = length(levels(WJ_Seurat_Obj$Development)))
  names(cellchat_list) <- levels(WJ_Seurat_Obj$Development)
  for(tp in levels(WJ_Seurat_Obj$Development)) {
    
    ### subset for the specific time point
    subset_obj <- subset(WJ_Seurat_Obj,
                         cells = rownames(WJ_Seurat_Obj@meta.data)[which(WJ_Seurat_Obj$Development == tp)])
    
    ### see how many cells in each cell group
    sapply(unique(subset_obj$computational_annotation), function(x) {
      return(length(which(subset_obj$computational_annotation == x)))
    })
    
    ### cell chat
    ### create a cell chat object
    cellchat <- createCellChat(object = subset_obj, group.by = "computational_annotation")
    
    ### set ligand-receptor interaction db
    ### CellChatDB <- CellChatDB.human
    CellChatDB <- CellChatDB.mouse
    
    ### use a subset of CellChatDB for cell-cell communication analysis
    ### CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
    ### use all CellChatDB for cell-cell communication analysis
    CellChatDB.use <- CellChatDB
    
    ### set the used database in the object
    cellchat@DB <- CellChatDB.use
    
    ### subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    
    ### Preprocessing
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    ### project gene expression data onto PPI network (optional)
    cellchat <- projectData(cellchat, PPI.mouse)
    
    ### Compute the communication probability and infer cellular communication network
    cellchat <- computeCommunProb(cellchat)
    ### Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 20)
    
    ### Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)
    
    ### Calculate the aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    
    ### save the cellchat result to the list
    cellchat_list[[tp]] <- cellchat
    
    gc()
    
  }
  
  ### make cellchat figures for each time point
  for(tp in levels(WJ_Seurat_Obj$Development)) {
    
    ### create output path for each time point
    outputDir2 <- paste0(outputDir, "/cellchat/", tp, "/")
    dir.create(path = outputDir2, recursive = TRUE)
    
    ### empty list
    plot_list[[tp]] <- list()
    
    ### We can also visualize the aggregated cell-cell communication network.
    ### For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
    groupSize <- as.numeric(table(cellchat_list[[tp]]@idents))
    png(filename = paste0(outputDir2, "/", tp, "_Cellchat_Circle_Interaction.png"), width = 3000, height = 2500, res = 350)
    netVisual_circle(cellchat_list[[tp]]@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()
    
    ### Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group.
    ### Here we also control the parameter edge.weight.max so that we can compare edge weights between different networks.
    mat <- cellchat_list[[tp]]@net$weight
    png(filename = paste0(outputDir2, "/", tp, "_Cellchat_Circle_Interaction_Each.png"), width = 4000, height = 3000, res = 300)
    par(mfrow = c(4,4), xpd=TRUE)
    for (i in 1:nrow(mat)) {
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[i, ] <- mat[i, ]
      netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,
                       edge.weight.max = max(mat), title.name = rownames(mat)[i],
                       vertex.label.cex = 0.4)
    }
    dev.off()
    
    ### One can visualize the inferred communication network of signaling pathways using netVisual_aggregate,
    ### and visualize the inferred communication networks of individual L-R pairs associated with that signaling pathway using netVisual_individual.
    print(cellchat_list[[tp]]@netP$pathways)
    pathways.show <- c("CXCL")
    png(filename = paste0(outputDir2, "/", tp, "_Cellchat_Circle_CXCL.png"), width = 3000, height = 2500, res = 350)
    netVisual_aggregate(cellchat_list[[tp]], signaling = pathways.show,  layout = "circle")
    dev.off()
    png(filename = paste0(outputDir2, "/", tp, "_Cellchat_Heatmap_CXCL.png"), width = 3000, height = 2500, res = 350)
    netVisual_heatmap(cellchat_list[[tp]], signaling = pathways.show, color.heatmap = "Reds")
    dev.off()
    
    ### Compute the contribution of each ligand-receptor pair to the overall signaling pathway and
    ### visualize cell-cell communication mediated by a single ligand-receptor pair
    png(filename = paste0(outputDir2, "/", tp, "_Cellchat_LR_Contribution_CXCL.png"), width = 3000, height = 2500, res = 350)
    netAnalysis_contribution(cellchat_list[[tp]], signaling = pathways.show)
    dev.off()
    
    ### show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use')
    ### to other cell groups (defined by 'targets.use')
    png(filename = paste0(outputDir2, "/", tp, "_Cellchat_Significant_LR.png"), width = 7000, height = 10000, res = 300)
    netVisual_bubble(cellchat_list[[tp]], sources.use = NULL, targets.use = NULL, remove.isolate = TRUE, thresh = 0.01)
    dev.off()
    
    ### Gene expression of the signaling pathway genes
    png(filename = paste0(outputDir2, "/", tp, "_Cellchat_Gene_Exp_CXCL.png"), width = 3000, height = 2500, res = 350)
    plotGeneExpression(cellchat, signaling = "CXCL")
    dev.off()
    
  }
  
  
}

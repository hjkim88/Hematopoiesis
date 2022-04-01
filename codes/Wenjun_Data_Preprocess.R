###
#   File name : Wenjun_Data_Preprocess.R
#   Author    : Hyunjin Kim
#   Date      : Feb 3, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Aggregate Cell Ranger output and preprocess the data
#   
#   Sample Annotation: 
#     TDH017            : P6 Heme
#     TDH018_2          : P14 Heme
#     TDH019_2          : P14 Stroma
#     TDH025_2          : P6 Stroma
#     TDH026_1350503    : P21 Heme (Trent 1st)
#     TDH026_1355923    : P21 Heme (Trent 2nd)
#     TDH027_1350504    : P21 Stroma (Trent 1st)
#     TDH027_1355924    : P21 Stroma (Trent 2nd)
#     Wenjun_Heme1_2    : P28 Heme
#     Wenjun_Heme2_2    : P21 Heme
#     Wenjun_Stromal1_2 : P28 Stroma
#     Wenjun_Stromal2_2 : P21 Stroma
#     TDH007            : P0 Heme
#     TDH008_1094045    : P0 Stroma (Trent 1st)
#     TDH008_1106032    : P0 Stroma (Trent 2nd)
#     TDH001            : ADULT Heme
#     TDH002            : ADULT Stroma
#
#   Instruction
#               1. Source("Wenjun_Data_Preprocess.R")
#               2. Run the function "wenjun_preprocess" - specify the input path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Wenjun_Data_Preprocess.R/Wenjun_Data_Preprocess.R")
#               > wenjun_preprocess(cellranger_result_dir="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Wenjun/Result/",
#                                   aggr_feature_bc_mat_dir="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Wenjun/Result/Wenjun_aggr3/outs/count/filtered_feature_bc_matrix/",
#                                   output_obj_path="./data/Wenjun_Seurat_Obj2.RDS")
###

wenjun_preprocess <- function(cellranger_result_dir="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Wenjun/Result/",
                              aggr_feature_bc_mat_dir="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Wenjun/Result/Wenjun_aggr3/outs/count/filtered_feature_bc_matrix/",
                              output_obj_path="./data/Wenjun_Seurat_Obj2.RDS") {
  
  ### load libraries
  if(!require(dplyr, quietly = TRUE)) {
    install.packages("dplyr")
    require(dplyr, quietly = TRUE)
  }
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(patchwork, quietly = TRUE)) {
    install.packages("patchwork")
    require(patchwork, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(grid, quietly = TRUE)) {
    install.packages("grid")
    require(grid, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  
  
  ### get separate folders in the cell ranger out dir
  sample_folders <- list.dirs(path = cellranger_result_dir, full.names = FALSE, recursive = FALSE)
  sample_folders <- sample_folders[-grep("aggr", sample_folders)]
  
  ### load metric files into R memory space
  metric_info <- vector("list", length = length(sample_folders))
  names(metric_info) <- sample_folders
  metric_info2 <- NULL
  for(sp in names(metric_info)) {
    writeLines(paste(sp))
    metric_info[[sp]] <- t(read.csv(file = paste0(cellranger_result_dir, sp, "/outs/metrics_summary.csv"),
                                    stringsAsFactors = FALSE, check.names = FALSE))
    if(is.null(metric_info2)) {
      metric_info2 <- metric_info[[sp]]
      colnames(metric_info2)[1] <- sp
    } else {
      metric_info2 <- merge(metric_info2, metric_info[[sp]], by = 0, all = TRUE)
      rownames(metric_info2) <- metric_info2[,1]
      metric_info2 <- metric_info2[,-1]
      colnames(metric_info2)[ncol(metric_info2)] <- sp
    }
  }
  
  ### draw bar graphs with some features
  plot_df <- data.frame(Sample=colnames(metric_info2),
                        Number_of_Reads=as.numeric(gsub(",", "", metric_info2["Number of Reads",])),
                        Number_of_Cells=as.numeric(gsub(",", "", metric_info2["Estimated Number of Cells",])),
                        Mean_Reads_per_Cell=as.numeric(gsub(",", "", metric_info2["Mean Reads per Cell",])),
                        Median_Genes_per_Cell=as.numeric(gsub(",", "", metric_info2["Median Genes per Cell",])),
                        Median_UMI_Counts_per_Cell=as.numeric(gsub(",", "", metric_info2["Median UMI Counts per Cell",])),
                        Total_Genes_Detected=as.numeric(gsub(",", "", metric_info2["Total Genes Detected",])),
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  p <- vector("list", length = ncol(plot_df)-1)
  names(p) <- colnames(plot_df)[-1]
  
  for(col in names(p)) {
    p[[col]] <- ggplot(data=plot_df, aes_string(x="Sample", y=col, label=col)) +
      geom_bar(stat = "identity", fill = "gray60") +
      ggtitle("") +
      geom_text(hjust=0.6, color="black", size=5) +
      coord_flip() +
      theme_classic(base_size = 30) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 36),
            axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
            axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 30),
            axis.ticks = element_blank())
  }
  
  g <- arrangeGrob(grobs = p,
                   nrow = 3,
                   ncol = 2)
  ggsave(file = paste0("./data/Wenjun_Data_QC_Results2.pdf"), g, width = 25, height = 22)
  
  
  ### load the aggregated Wenjun dataset generated by cellranger aggr
  wenjun.data <- Read10X(data.dir = aggr_feature_bc_mat_dir)
  
  ### create a Seurat object
  ### min.cells : include genes detected in at least this many cells
  ### min.features : include cells where at least this many features are detected
  Wenjun_Seurat_Obj <- CreateSeuratObject(counts = wenjun.data,
                                          project = "Wenjun_HSPC",
                                          min.cells = 3,
                                          min.features = 200)
  
  ### dissect the combined barcodes
  Wenjun_Seurat_Obj$orig.barcode <- sapply(rownames(Wenjun_Seurat_Obj@meta.data), function(x) {
    return(strsplit(x, split = "-", fixed = TRUE)[[1]][1])
  })
  Wenjun_Seurat_Obj$lib.ident <- sapply(rownames(Wenjun_Seurat_Obj@meta.data), function(x) {
    return(strsplit(x, split = "-", fixed = TRUE)[[1]][2])
  })
  
  ### jacard index
  cal_jaccard <- function(x, y) {
    return(length(intersect(x, y)) / length(union(x, y)))
  }
  
  ### annotate the cells which are which
  ### the barcodes don't match exactly because when creating seurat, there was a filtering process
  Wenjun_Seurat_Obj$library <- NA
  for(sample_name in sample_folders) {
    
    temp_barcodes <- as.character(read.csv(gzfile(paste0(cellranger_result_dir, sample_name, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")),
                                           header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)[,1])
    temp_barcodes <- sapply(temp_barcodes, function(x) {
      return(strsplit(x, split = "-", fixed = TRUE)[[1]][1])
    }, USE.NAMES = FALSE)
    
    jaccard_list <- rep(0, length(unique(Wenjun_Seurat_Obj$lib.ident)))
    names(jaccard_list) <- unique(Wenjun_Seurat_Obj$lib.ident)
    for(lib in unique(Wenjun_Seurat_Obj$lib.ident)) {
      lib_barcodes <- Wenjun_Seurat_Obj$orig.barcode[which(Wenjun_Seurat_Obj$lib.ident == lib)]
      
      jaccard_list[lib] <- cal_jaccard(temp_barcodes, lib_barcodes)
    }
    
    if(length(which(jaccard_list > 0.8)) > 1) {
      writeLines(paste("ERROR: There are multiple samples matched with this barcodes:", sample_name))
    }
    
    Wenjun_Seurat_Obj$library[which(Wenjun_Seurat_Obj$lib.ident == names(jaccard_list)[which(jaccard_list > 0.8)])] <- sample_name
    
  }
  
  ### remove TDH026_1350503, TDH026_1355923, TDH027_1350504, TDH027_1355924
  ### since they are low quality data
  Wenjun_Seurat_Obj <- subset(Wenjun_Seurat_Obj,
                              cells = rownames(Wenjun_Seurat_Obj@meta.data)[-which(Wenjun_Seurat_Obj$library %in% c("TDH026_1350503",
                                                                                                                    "TDH026_1355923",
                                                                                                                    "TDH027_1350504",
                                                                                                                    "TDH027_1355924"))])
  
  ### annotate some more
  Wenjun_Seurat_Obj$Development <- ""
  Wenjun_Seurat_Obj$Development[which(Wenjun_Seurat_Obj$library %in% c("TDH017", "TDH025_2"))] <- "P6"
  Wenjun_Seurat_Obj$Development[which(Wenjun_Seurat_Obj$library %in% c("TDH018_2", "TDH019_2"))] <- "P14"
  Wenjun_Seurat_Obj$Development[which(Wenjun_Seurat_Obj$library %in% c("Wenjun_Heme2_2", "Wenjun_Stromal2_2"))] <- "P21"
  Wenjun_Seurat_Obj$Development[which(Wenjun_Seurat_Obj$library %in% c("Wenjun_Heme1_2", "Wenjun_Stromal1_2"))] <- "P28"
  Wenjun_Seurat_Obj$Development[which(Wenjun_Seurat_Obj$library %in% c("TDH007", "TDH008_1094045", "TDH008_1106032"))] <- "P0"
  Wenjun_Seurat_Obj$Development[which(Wenjun_Seurat_Obj$library %in% c("TDH001", "TDH002"))] <- "ADULT"
  
  Wenjun_Seurat_Obj$Cell1 <- ""
  Wenjun_Seurat_Obj$Cell1[which(Wenjun_Seurat_Obj$library %in% c("TDH017", "TDH018_2", "Wenjun_Heme2_2", "Wenjun_Heme1_2", "TDH007", "TDH001"))] <- "Heme"
  Wenjun_Seurat_Obj$Cell1[which(Wenjun_Seurat_Obj$library %in% c("TDH019_2", "TDH025_2", "Wenjun_Stromal2_2", "Wenjun_Stromal1_2", "TDH008_1094045", "TDH008_1106032", "TDH002"))] <- "Stroma"
  
  Wenjun_Seurat_Obj$Cell2 <- paste0(Wenjun_Seurat_Obj$Development, "_", Wenjun_Seurat_Obj$Cell1)
  
  ### MT percentage (Human: ^MT-, Mouse: ^mt- or ^Mt-)
  Wenjun_Seurat_Obj[["percent.mt"]] <- PercentageFeatureSet(Wenjun_Seurat_Obj, pattern = "^mt-")
  
  ### Visualize QC metrics as a violin plot
  VlnPlot(Wenjun_Seurat_Obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot(density(Wenjun_Seurat_Obj@meta.data$nFeature_RNA))
  axis(1, seq(0,8000,500), tck = 1, col = "gray")
  plot(density(Wenjun_Seurat_Obj@meta.data$nCount_RNA))
  axis(1, seq(0,100000,5000), tck = 1, col = "gray")
  plot(density(Wenjun_Seurat_Obj@meta.data$percent.mt))
  FeatureScatter(Wenjun_Seurat_Obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
  FeatureScatter(Wenjun_Seurat_Obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
  
  ### Filter cells that have unique feature counts over 7000 or less than 1000
  ### Filter cells that have RNA mocule # > 50000 or RNA mocule # < 3000
  ### Filter cells that have > 10% mitochondrial counts
  Wenjun_Seurat_Obj <- subset(Wenjun_Seurat_Obj, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & nCount_RNA > 3000 & nCount_RNA < 50000 & percent.mt < 10)
  
  ### normalization
  Wenjun_Seurat_Obj <- NormalizeData(Wenjun_Seurat_Obj,
                                     normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### Basic function to convert human to mouse gene names
  convertHumanGeneList <- function(x){
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
    humanx <- unique(genesV2[, 2])
    
    return(humanx)
  }
  
  ### Cell cycle score (will be used later for regression out)
  Wenjun_Seurat_Obj <- CellCycleScoring(object = Wenjun_Seurat_Obj,
                                        g2m.features = convertHumanGeneList(cc.genes$g2m.genes),
                                        s.features = convertHumanGeneList(cc.genes$s.genes))
  
  ### find variable genes
  Wenjun_Seurat_Obj <- FindVariableFeatures(Wenjun_Seurat_Obj,
                                            selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  Wenjun_Seurat_Obj <- ScaleData(Wenjun_Seurat_Obj,
                                 vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
  
  ### PCA
  Wenjun_Seurat_Obj <- RunPCA(Wenjun_Seurat_Obj,
                              features = VariableFeatures(object = Wenjun_Seurat_Obj), npcs = 15)
  ElbowPlot(Wenjun_Seurat_Obj, ndims = 15, reduction = "pca")
  
  ### UMAP
  Wenjun_Seurat_Obj <- RunUMAP(Wenjun_Seurat_Obj, dims = 1:15)
  
  ### clustering
  Wenjun_Seurat_Obj <- FindNeighbors(Wenjun_Seurat_Obj, dims = 1:15)
  Wenjun_Seurat_Obj <- FindClusters(Wenjun_Seurat_Obj, resolution = 0.1)
  
  ### save the Seurat object as RDS file
  saveRDS(Wenjun_Seurat_Obj, file = output_obj_path)
  
  ### factorize some columns
  Wenjun_Seurat_Obj$Development <- factor(Wenjun_Seurat_Obj$Development, levels = c("P0", "P6", "P14", "P21", "P28", "ADULT"))
  Wenjun_Seurat_Obj$Cell2 <- factor(Wenjun_Seurat_Obj$Cell2,
                                    levels = c("P0_Heme", "P0_Stroma",
                                               "P6_Heme", "P6_Stroma",
                                               "P14_Heme", "P14_Stroma",
                                               "P21_Heme", "P21_Stroma",
                                               "P28_Heme", "P28_Stroma",
                                               "ADULT_Heme", "ADULT_Stroma"))
  
  ### draw a UMAP
  p <- list()
  p[[1]] <- DimPlot(object = Wenjun_Seurat_Obj, reduction = "umap", raster = TRUE,
               group.by = "Development",
               pt.size = 0.5) +
    ggtitle(paste0("")) +
    labs(color = "Development") +
    scale_color_brewer(palette="Set1") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 48) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 36, color = "black", face = "bold"),
          axis.title = element_text(size = 36, color = "black", face = "bold"),
          axis.text = element_text(size = 36, color = "black", face = "bold"),
          legend.title = element_text(size = 24, color = "black", face = "bold"),
          legend.text = element_text(size = 24, color = "black", face = "bold"),
          axis.ticks = element_blank())
  p[[2]] <- DimPlot(object = Wenjun_Seurat_Obj, reduction = "umap", raster = TRUE,
                    group.by = "Cell1",
                    pt.size = 0.5) +
    ggtitle(paste0("")) +
    labs(color = "Heme/Stroma") +
    scale_color_brewer(palette="Set1") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 48) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 36, color = "black", face = "bold"),
          axis.title = element_text(size = 36, color = "black", face = "bold"),
          axis.text = element_text(size = 36, color = "black", face = "bold"),
          legend.title = element_text(size = 24, color = "black", face = "bold"),
          legend.text = element_text(size = 24, color = "black", face = "bold"),
          axis.ticks = element_blank())
  p[[3]] <- DimPlot(object = Wenjun_Seurat_Obj, reduction = "umap", raster = TRUE,
                    group.by = "Cell2",
                    pt.size = 0.5) +
    ggtitle(paste0("")) +
    labs(color = "Cell Type") +
    scale_color_brewer(palette="Set1") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 48) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 36, color = "black", face = "bold"),
          axis.title = element_text(size = 36, color = "black", face = "bold"),
          axis.text = element_text(size = 36, color = "black", face = "bold"),
          legend.title = element_text(size = 24, color = "black", face = "bold"),
          legend.text = element_text(size = 24, color = "black", face = "bold"),
          axis.ticks = element_blank())
  p[[4]] <- DimPlot(object = Wenjun_Seurat_Obj, reduction = "umap", raster = TRUE,
                    group.by = "seurat_clusters",
                    pt.size = 0.5) +
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
  
  ### print out the UMAPs
  g <- arrangeGrob(grobs = p,
                   nrow = 2,
                   ncol = 2,
                   top = "")
  ggsave(paste0("./data/Wenjun_Total_UMAP.png"), plot = g, width = 30, height = 16, dpi = 350)
  
  ggsave(paste0("./data/Wenjun_UMAP_Development.png"), plot = p[[1]], width = 20, height = 10, dpi = 350)
  ggsave(paste0("./data/Wenjun_UMAP_Cell_Type1.png"), plot = p[[2]], width = 20, height = 10, dpi = 350)
  ggsave(paste0("./data/Wenjun_UMAP_Cell_Type2.png"), plot = p[[3]], width = 20, height = 10, dpi = 350)
  ggsave(paste0("./data/Wenjun_UMAP_Clusters.png"), plot = p[[4]], width = 20, height = 10, dpi = 350)
  
  ### print some info for Shannon
  sapply(levels(Wenjun_Seurat_Obj$Development), function(x) {
    return(length(intersect(which(Wenjun_Seurat_Obj$Development == x),
                            which(Wenjun_Seurat_Obj$Cell1 == "Heme"))))
  })
  sapply(levels(Wenjun_Seurat_Obj$Development), function(x) {
    return(length(intersect(which(Wenjun_Seurat_Obj$Development == x),
                            which(Wenjun_Seurat_Obj$Cell1 == "Stroma"))))
  })
  
  
  
}
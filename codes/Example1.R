### THIS IS AN EXAMPLE CODE FOR PC1 GENE - HEATMAP & PATHWAY ANALYSIS

### set parameters
Robj_path <- "./data/Combined_Seurat_Obj.RDATA"
outputDir <- "./results/Pseudotime/"

### load the Seurat object and save the object name
tmp_env <- new.env()
load(paste0(Robj_path), tmp_env)
obj_name <- ls(tmp_env)
assign("Combined_Seurat_Obj", get(obj_name, envir = tmp_env))
rm(tmp_env)
gc()

### run PCA
Combined_Seurat_Obj <- FindVariableFeatures(Combined_Seurat_Obj)
Combined_Seurat_Obj <- ScaleData(Combined_Seurat_Obj)
Combined_Seurat_Obj <- RunPCA(Combined_Seurat_Obj, npcs = 10)

### draw a PCA
DimPlot(Combined_Seurat_Obj, reduction = "pca", group.by = "Tissue", pt.size = 1.5) +
  labs(title = paste0("PCA_Combined_Tissue"))

### check whether the orders are the same
print(identical(names(Idents(object = Combined_Seurat_Obj)), rownames(Combined_Seurat_Obj@meta.data)))

### determine necessary variables
time_points <- c("E16", "E18", "P0", "ADULT")
HSPC_populations <- c("LTHSC", "STHSC", "MPP2", "MPP3", "MPP4")
possible_names <- as.vector(sapply(time_points, function(x) paste0(x, HSPC_populations)))
possible_names_mat <- data.frame(Names=possible_names,
                                 Time=as.vector(sapply(time_points, function(x) rep(x, length(HSPC_populations)))),
                                 HSPC=rep(HSPC_populations, length(time_points)),
                                 stringsAsFactors = FALSE, check.names = FALSE)
row.names(possible_names_mat) <- possible_names

#
### separate the combined into different cell groups
#

### set the ident of the object with the HSPC type
Combined_Seurat_Obj <- SetIdent(object = Combined_Seurat_Obj,
                                cells = rownames(Combined_Seurat_Obj@meta.data),
                                value = Combined_Seurat_Obj@meta.data$HSPC)

### check whether the orders are the same
print(identical(names(Idents(object = Combined_Seurat_Obj)), rownames(Combined_Seurat_Obj@meta.data)))

### rownames in the meta.data should be in the same order as colnames in the counts
Combined_Seurat_Obj@meta.data <- Combined_Seurat_Obj@meta.data[colnames(Combined_Seurat_Obj@assays$RNA@counts),]

### active assay = "RNA"
Combined_Seurat_Obj@active.assay <- "RNA"

### choose cell type
type <- "MPP2"

### new output directory
outputDir2 <- paste0(outputDir, type, "/")
dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)

### split the Seurat obj based on HSPC info
Seurat_Obj <- subset(Combined_Seurat_Obj, idents=type)

### order the meta data by developmental time
Seurat_Obj@meta.data <- Seurat_Obj@meta.data[order(factor(Seurat_Obj@meta.data$Development,
                                                          levels = time_points)),]

### rownames in the meta.data should be in the same order as colnames in the counts
Seurat_Obj@assays$RNA@counts <- Seurat_Obj@assays$RNA@counts[,rownames(Seurat_Obj@meta.data)]

### run PCA
Seurat_Obj <- RunPCA(Seurat_Obj, npcs = 10)
pca_map <- Embeddings(Seurat_Obj, reduction = "pca")[rownames(Seurat_Obj@meta.data),1:10]


### find feature contributions of the PC1 
pca_cos2 <- Seurat_Obj@reductions$pca@feature.loadings * Seurat_Obj@reductions$pca@feature.loadings
pca_contb <- pca_cos2
for(i in 1:ncol(pca_contb)) {
  s <- sum(pca_cos2[,i])
  for(j in 1:nrow(pca_contb)) {
    pca_contb[j,i] <- pca_cos2[j,i] * 100 / s
  }
}
pca_contb <- pca_contb[order(-pca_contb[,"PC_1"]),]

### get genes that contributed to the PC1 the most
contb_threshold <- 0.1
important_genes <- rownames(pca_contb)[which(pca_contb[,"PC_1"] > contb_threshold)]

#
### a heatmap with the genes - 322 genes x 428 cells
#

### get a matrix for the heatmap
heatmap_mat <- data.frame(Seurat_Obj@assays$RNA@counts[important_genes,], check.names = FALSE)

### A function for scaling for heatmap
scale_h <- function(data, type, na.rm=TRUE) {
  
  if(type == "row") {
    scaled <- t(scale(t(data)))
  } else if(type == "col") {
    scaled <- scale(data)
  } else {
    stop("Type is required: row or col")
  }
  
  if(na.rm == TRUE && (length(which(is.na(scaled))) > 0))  {
    scaled <- scaled[-unique(which(is.na(scaled), arr.ind = TRUE)[,1]),]
  }
  
  return(scaled)
}

### scale the data
heatmap_mat_scaled <- scale_h(heatmap_mat, type = "row")

### see the distribution of the gene expressions
# plot(density(heatmap_mat_scaled))

### because there are some outliers in positive values
### we set the maximum as abs(minimum)
heatmap_mat_scaled[which(heatmap_mat_scaled > abs(min(heatmap_mat_scaled)))] <- abs(min(heatmap_mat_scaled))

### set colside colors
uniqueV <- unique(Seurat_Obj@meta.data$Development)
colors <- colorRampPalette(brewer.pal(9,"Blues"))(length(uniqueV))
names(colors) <- uniqueV

### hierarchical clustering functions
dist.spear <- function(x) as.dist(1-cor(t(x), method = "spearman"))
hclust.ave <- function(x) hclust(x, method="average")

### heatmap
png(paste0(outputDir2, "PC1_Genes_Heatmap.png"), width = 2000, height = 1000)
par(oma=c(0,0,2,6))
heatmap.3(as.matrix(heatmap_mat_scaled), main = paste0("PC1_Genes_Heatmap_(",
                                                       nrow(heatmap_mat_scaled), " Genes x ",
                                                       ncol(heatmap_mat_scaled), " Cells)"),
          xlab = "", ylab = "", col=greenred(300),
          scale="none", key=T, keysize=0.8, density.info="density",
          dendrogram = "none", trace = "none",
          labRow = rownames(heatmap_mat_scaled), labCol = FALSE,
          Rowv = TRUE, Colv = FALSE,
          distfun=dist.spear, hclustfun=hclust.ave,
          ColSideColors = cbind(colors[as.character(Seurat_Obj@meta.data$Development)]),
          cexRow = 2, cexCol = 2, na.rm = TRUE)
legend("left", inset = 0, xpd = TRUE, title = "Time", legend = names(colors), fill = colors, cex = 2, box.lty = 0)
dev.off()

### pathway analysis with the important genes of the PC1
contb_threshold <- 0.1
important_genes <- rownames(pca_contb)[which(pca_contb[,"PC_1"] > contb_threshold)]
pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Mm.eg.db,
                                                          important_genes,
                                                          "ENTREZID", "SYMBOL"),
                                        org = "mouse", database = "GO",
                                        title = paste0("Pathway_Results_", length(important_genes), "_PC1_Genes_", contb_threshold),
                                        displayNum = 50, imgPrint = TRUE,
                                        dir = paste0(outputDir2))
pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Mm.eg.db,
                                                            important_genes,
                                                            "ENTREZID", "SYMBOL"),
                                          org = "mouse", database = "KEGG",
                                          title = paste0("Pathway_Results_", length(important_genes), "_PC1_Genes_", contb_threshold),
                                          displayNum = 50, imgPrint = TRUE,
                                          dir = paste0(outputDir2))
write.xlsx2(pathway_result_GO, file = paste0(outputDir2, "GO_pathway_results_", length(important_genes), "_PC1_Genes_", contb_threshold, ".xlsx"),
            row.names = FALSE, sheetName = paste0("GO_Results"))
write.xlsx2(pathway_result_KEGG, file = paste0(outputDir2, "KEGG_pathway_results_", length(important_genes), "_PC1_Genes_", contb_threshold, ".xlsx"),
            row.names = FALSE, sheetName = paste0("KEGG_Results"))

### Analyses for Trent - 04/26/2021

### monocle plot_complex_cell_trajectory modify

### set paths
Robj1_path="./data/Combined_Seurat_Obj.RDATA"
Robj2_path="./data/Combined_Seurat_Obj.RDS"
outputDir="./results/Apr_2021/"

dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

### load libraries
if(!require(Seurat, quietly = TRUE)) {
  install.packages("Seurat")
  require(Seurat, quietly = TRUE)
}
if(!require(ggplot2, quietly = TRUE)) {
  install.packages("ggplot2")
  require(ggplot2, quietly = TRUE)
}
if(!require(ggrepel, quietly = TRUE)) {
  install.packages("ggrepel")
  require(ggrepel, quietly = TRUE)
}
if(!require(slingshot, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("slingshot")
  require(slingshot, quietly = TRUE)
}
if(!require(scales, quietly = TRUE)) {
  install.packages("scales")
  library(scales, quietly = TRUE)
}
if(!require(RColorBrewer, quietly = TRUE)) {
  install.packages("RColorBrewer")
  require(RColorBrewer, quietly = TRUE)
}
if(!require(gplots, quietly = TRUE)) {
  install.packages("gplots")
  library(gplots, quietly = TRUE)
}
if(!require(msigdbr, quietly = TRUE)) {
  install.packages("msigdbr")
  library(msigdbr, quietly = TRUE)
}
if(!require(xlsx, quietly = TRUE)) {
  install.packages("xlsx")
  require(xlsx, quietly = TRUE)
}
if(!require(org.Mm.eg.db, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("org.Mm.eg.db")
  require(org.Mm.eg.db, quietly = TRUE)
}
if(!require(ggupset, quietly = TRUE)) {
  install.packages("ggupset")
  require(ggupset, quietly = TRUE)
}
if(!require(ggnewscale, quietly = TRUE)) {
  install.packages("ggnewscale")
  require(ggnewscale, quietly = TRUE)
}
if(!require(ggthemes, quietly = TRUE)) {
  install.packages("ggthemes")
  require(ggthemes, quietly = TRUE)
}
if(!require(msigdbr, quietly = TRUE)) {
  install.packages("msigdbr")
  library(msigdbr, quietly = TRUE)
}
if(!require(DOSE, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DOSE")
  require(DOSE, quietly = TRUE)
}
if(!require(enrichplot, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("enrichplot")
  require(enrichplot, quietly = TRUE)
}
if(!require(monocle, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("monocle")
  require(monocle, quietly = TRUE)
}


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
              theme_classic(base_size = 50) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
              scale_x_discrete(limits = rev(description)) +
              guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
              ggtitle(paste0("KEGG ", title)) +
              theme(axis.text = element_text(size = 50))
            
            png(paste0(dir, "kegg_", title, "_CB.png"), width = 2200, height = 1000)
            print(p[[1]])
            dev.off()
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
              theme_classic(base_size = 50) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
              scale_x_discrete(limits = rev(description)) +
              guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
              ggtitle(paste0("GO ", title)) +
              theme(axis.text = element_text(size = 50))
            
            png(paste0(dir, "go_", title, "_CB.png"), width = 2200, height = 1000)
            print(p[[2]])
            dev.off()
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

#
### GSEA with the important genes of the PC1
#

#'****************************************************************************************
#' Gene Set Enrichment Analysis function
#' 
#' It receives gene list (character vector) and signature profiles (named numeric vector)
#' as inputs, performs GSEA and returns a table of GSEA result table and draws
#' a GSEA plot. It is basically a statistical significance test to check how the
#' given given genes are biased on the signature profiles.
#' 
#' Whether there are multiple gene sets or multiple signatures,
#' multiple testing (FDR computation) is performed.
#' But if the input gene set and the input signature are both lists with multiple
#' items (The length of the two are both more than 1) then we return an error message.
#' 
#' The plot file names will be determined by names(gene_list) or names(signature)
#' If length(gene_list) > 1, then names(gene_list) will be used and
#' if length(signature) > 1, then names(signature) will be used as file names.
#' If there is no list names, then file names will be "GSEA_Plot_i.png".
#' Here, i indicates that the plot is from the i-th row of the GSEA result table.
#' 
#' * Some plot drawing codes were from Rtoolbox/R/ReplotGSEA.R written by Thomas Kuilman. 
#'****************************************************************************************
#' @title	run_gsea
#' 
#' @param gene_list   A list of character vectors containing gene names to be tested
#' @param signature   A list of named numeric vectors of signature values for GSEA. The gene_list
#'                    should be included in the names(signature)
#' @param printPlot   If TRUE, it also generates GSEA plot of the results
#'                    (Default = FALSE)
#' @param fdr_cutoff  When printing GSEA plots, print them with the FDR < fdr_cutoff only
#'                    (Default = 0.05)
#' @param printPath   When printing GSEA plots, print them in the designated path
#'                    (Default = "./")
#' @param width       The width of the plot file
#'                    (Default = 2000)
#' @param height      The height of the plot file
#'                    (Default = 1200)
#' @param res         The resolution of the plot file
#'                    (Default = 130)
#' 
#' @return 	          It tests bias of the "gene_list" on the "signature" range and
#'                    returns a table including p-values and FDRs (adjusted p-values)
#'                    If fdr_cutoff == TRUE, it also generates a GSEA plot with the result
#' 
run_gsea <- function(gene_list,
                     signature,
                     printPlot = FALSE,
                     fdr_cutoff = 0.05,
                     width = 2000,
                     height = 1200,
                     res = 130,
                     printPath = "./") {
  
  ### load required libraries
  if(!require("fgsea", quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("fgsea")
    require("fgsea", quietly = TRUE)
  }
  if(!require("limma", quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("limma")
    require("limma", quietly = TRUE)
  } 
  if(!require("checkmate", quietly = TRUE)) {
    install.packages("checkmate")
    require("checkmate", quietly = TRUE)
  }
  
  ### argument checking
  assertList(gene_list)
  assertList(signature)
  assertLogical(printPlot)
  assertNumeric(fdr_cutoff)
  assertIntegerish(width)
  assertIntegerish(height)
  assertIntegerish(res)
  assertString(printPath)
  if(length(gene_list) > 1 && length(signature) > 1) {
    stop("ERROR: \"gene_list\" and \"signature\" cannot be both \"list\"")
  }
  
  ### set random seed
  set.seed(1234)
  
  ### run GSEA
  ### if there are more than one signatures
  if(length(signature) > 1) {
    ### combine GSEA results of every signature inputs
    for(i in 1:length(signature)) {
      temp <- data.frame(fgsea(pathways = gene_list, stats = signature[[i]], nperm = 1000))
      if(i == 1) {
        gsea_result <- temp
      } else {
        gsea_result <- rbind(gsea_result, temp)
      }
    }
    
    ### compute FDRs
    corrected_gsea_result <- gsea_result[order(gsea_result$pval),]
    corrected_gsea_result$padj <- p.adjust(corrected_gsea_result$pval, method = "BH")
    gsea_result <- corrected_gsea_result[rownames(gsea_result),]
  }
  ### if there are more than one gene sets
  else {
    gsea_result <- data.frame(fgsea(pathways = gene_list, stats = signature[[1]], nperm = 1000))
  }
  
  ### print GSEA plot
  sIdx <- which(gsea_result$padj < fdr_cutoff)
  if(printPlot && length(sIdx) > 0) {
    for(i in sIdx) {
      ### get required values ready
      if(length(signature) > 1) {
        geneset <- gene_list[[1]]
        stats <- signature[[i]]
        stats <- stats[order(-stats)]
        fileName <- names(signature)[i]
      } else {
        geneset <- gene_list[[i]]
        stats <- signature[[1]]
        stats <- stats[order(-stats)]
        fileName <- names(gene_list)[i]
      }
      if(is.null(fileName)) {
        fileName <- paste0("GSEA_Plot_", i)
      }
      stats <- stats[!is.na(stats)]
      gsea.hit.indices <- which(names(stats) %in% geneset)
      es.temp <- calcGseaStat(stats, gsea.hit.indices, returnAllExtremes = TRUE)
      if(es.temp$res >= 0) {
        gsea.es.profile <- es.temp$tops
      } else {
        gsea.es.profile <- es.temp$bottoms
      }
      enrichment.score.range <- c(min(gsea.es.profile), max(gsea.es.profile))
      metric.range <- c(min(stats), max(stats))
      gsea.p.value <- round(gsea_result$pval[i] ,5)
      gsea.fdr <- round(gsea_result$padj[i] ,5)
      gsea.enrichment.score <- round(gsea_result$ES[i], 5)
      gsea.normalized.enrichment.score <- round(gsea_result$NES[i], 5)
      
      ### print GSEA result plot
      png(paste0(printPath, fileName, ".png"), width = width, height = height, res = res)
      
      ### set layout
      layout.show(layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.5, 0.2, 2)))
      
      ### draw the GSEA plot
      par(mar = c(0, 5, 2, 2))
      plot(c(1, gsea.hit.indices, length(stats)),
           c(0, gsea.es.profile, 0), type = "l", col = "red", lwd = 1.5, xaxt = "n",
           xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
           ylim = enrichment.score.range,
           main = list(fileName, font = 1, cex = 1),
           panel.first = {
             abline(h = seq(round(enrichment.score.range[1], digits = 1),
                            enrichment.score.range[2], 0.1),
                    col = "gray95", lty = 2)
             abline(h = 0, col = "gray50", lty = 2)
           }
      )
      
      ### add informative text to the GSEA plot
      plot.coordinates <- par("usr")
      if(es.temp$res < 0) {
        text(length(stats) * 0.01, plot.coordinates[3] * 0.98,
             paste("P-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
                   gsea.enrichment.score, "\nNormalized ES:",
                   gsea.normalized.enrichment.score), adj = c(0, 0))
      } else {
        text(length(stats) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
             paste("P-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
                   gsea.enrichment.score, "\nNormalized ES:",
                   gsea.normalized.enrichment.score), adj = c(1, 1))
      }
      
      ### draw hit indices
      par(mar = c(0, 5, 0, 2))
      plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
           ylab = "", xlim = c(1, length(stats)))
      abline(v = gsea.hit.indices, lwd = 0.75)
      
      ### create color palette for the heatmap
      par(mar = c(0, 5, 0, 2))
      rank.colors <- stats - metric.range[1]
      rank.colors <- rank.colors / (metric.range[2] - metric.range[1])
      rank.colors <- ceiling(rank.colors * 511 + 1)
      rank.colors <- colorRampPalette(c("blue", "white", "red"))(512)[rank.colors]
      
      ### draw the heatmap
      rank.colors <- rle(rank.colors)
      barplot(matrix(rank.colors$lengths), col = rank.colors$values,
              border = NA, horiz = TRUE, xaxt = "n", xlim = c(1, length(stats)))
      box()
      text(length(stats) / 2, 0.7,
           labels = "Signature")
      text(length(stats) * 0.01, 0.7, "Largest", adj = c(0, NA))
      text(length(stats) * 0.99, 0.7, "Smallest", adj = c(1, NA))
      
      ### draw signature values
      par(mar = c(5, 5, 0, 2))
      rank.metric <- rle(round(stats, digits = 2))
      plot(stats, type = "n", xaxs = "i",
           xlab = "Rank in ordered gene list", xlim = c(0, length(stats)),
           ylim = metric.range, yaxs = "i",
           ylab = "Signature values",
           panel.first = abline(h = seq(metric.range[1] / 2,
                                        metric.range[2] - metric.range[1] / 4,
                                        metric.range[2] / 2), col = "gray95", lty = 2))
      
      barplot(rank.metric$values, col = "lightgrey", lwd = 0.1,
              xlim = c(0, length(stats)), ylim = c(-1, 1),
              width = rank.metric$lengths, border = NA,
              space = 0, add = TRUE, xaxt = "n")
      box()
      
      ### print out the file
      dev.off()
    }
  }
  
  return(gsea_result)
  
}

### a function to get important & garbage genes and pathway analysis with those genes, GSEA with PC contribution,
### pathway analysis with DE genes from comparison, GSEA with DE genes from comparison.
### Seurat_object: The seurat object to be analysed
### target ident: the ident that will be analyzed, if NULL, just use the given object without split
### target_col: a column name of the meta data of the seurat object that will be used in
###             the pseudotime analysis and in creating heatmaps. Usually a time-associated column
### target_col_factor_level: a factor level of the 'target_col'
### one_size_only: if TRUE, only up & down-regulated genes will be used
### PC: the principle component that will be analysed, default: PC1
### PC_Val: the value of the given pc that will be used as a cutoff
###         the comparison will be made based on this value
### important_thresh: a contribution cut-off that will tell which genes contributed to the PC the most
### garbage_thresh: a contribution cut-off that will tell which genes contributed to the PC the least
### result_dir: the directory that all the results will be stored 
multiple_analyses_in_one <- function(Seurat_Object,
                                     target_ident=NULL,
                                     target_col,
                                     target_col_factor_level,
                                     one_side_only=FALSE,
                                     species=c("human", "mouse"),
                                     PC="PC_1",
                                     PC_Val=NULL,
                                     important_thresh=0.1,
                                     garbage_thresh=1e-04,
                                     result_dir="./") {
  
  ### check a struture
  if(!identical(names(Idents(object = Combined_Seurat_Obj)), rownames(Combined_Seurat_Obj@meta.data))) {
    stop("ERROR: Order of Idents of the given Seurat object does not match to that of the meta data")
  }
  
  ### create the result directory if it does not exist
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
  
  ### split the Seurat obj based on the given info
  if(is.null(target_ident)) {
    local_Seurat_Obj <- Seurat_Object
  } else {
    local_Seurat_Obj <- subset(Seurat_Object, idents=target_ident)
  }
  
  ### order the meta data by developmental time
  local_Seurat_Obj@meta.data <- local_Seurat_Obj@meta.data[order(factor(local_Seurat_Obj@meta.data[,target_col],
                                                                        levels = target_col_factor_level)),]
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  local_Seurat_Obj@assays$RNA@counts <- local_Seurat_Obj@assays$RNA@counts[,rownames(local_Seurat_Obj@meta.data)]
  
  ### run PCA
  local_Seurat_Obj <- RunPCA(local_Seurat_Obj, npcs = 10)
  pca_map <- Embeddings(local_Seurat_Obj, reduction = "pca")[rownames(local_Seurat_Obj@meta.data),1:10]
  
  ### find feature contributions of the given PC
  pca_cos2 <- local_Seurat_Obj@reductions$pca@feature.loadings * local_Seurat_Obj@reductions$pca@feature.loadings
  pca_contb <- pca_cos2
  for(i in 1:ncol(pca_contb)) {
    s <- sum(pca_cos2[,i])
    for(j in 1:nrow(pca_contb)) {
      pca_contb[j,i] <- pca_cos2[j,i] * 100 / s
    }
  }
  pca_contb <- pca_contb[order(-pca_contb[,PC]),]
  
  ### write out the PC contributions
  write.xlsx2(data.frame(Gene_Symbol=rownames(pca_contb),
                         pca_contb,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(result_dir, PC, "_Genes_Contributions.xlsx"),
              sheetName = paste0(PC, "_Genes_Contributions"),
              row.names = FALSE)
  
  ### get genes that contributed to the PC1 the most
  important_genes <- rownames(pca_contb)[which(pca_contb[,PC] > important_thresh)]
  
  ### get genes that contributed to the PC1 the least
  garbage_genes <- rownames(pca_contb)[which(pca_contb[,PC] < garbage_thresh)]
  
  ### pathway analysis with the important genes of the given PC
  if(species[1] == "human") {
    ### get entrez ids for the genes
    important_entrez_ids <- mapIds(org.Hs.eg.db,
                                   important_genes,
                                   "ENTREZID", "SYMBOL")
    important_entrez_ids <- important_entrez_ids[!is.na(important_entrez_ids)]
    
    pathway_result_GO <- pathwayAnalysis_CP(geneList = important_entrez_ids,
                                            org = species[1], database = "GO",
                                            title = paste0(PC, "_Pathway_Results"),
                                            displayNum = 10, imgPrint = TRUE,
                                            dir = paste0(result_dir))
    pathway_result_KEGG <- pathwayAnalysis_CP(geneList = important_entrez_ids,
                                              org = species[1], database = "KEGG",
                                              title = paste0(PC, "_Pathway_Results"),
                                              displayNum = 10, imgPrint = TRUE,
                                              dir = paste0(result_dir))
    
    ### enrichment object
    eo <- enrichGO(gene = important_entrez_ids, OrgDb = 'org.Hs.eg.db', readable = T, ont = "BP", pvalueCutoff = 0.05)
    eox <- setReadable(eo, 'org.Hs.eg.db', 'ENTREZID')
  } else if(species[1] == "mouse") {
    ### get entrez ids for the genes
    important_entrez_ids <- mapIds(org.Mm.eg.db,
                                   important_genes,
                                   "ENTREZID", "SYMBOL")
    important_entrez_ids <- important_entrez_ids[!is.na(important_entrez_ids)]
    
    pathway_result_GO <- pathwayAnalysis_CP(geneList = important_entrez_ids,
                                            org = species[1], database = "GO",
                                            title = paste0(PC, "_Pathway_Results"),
                                            displayNum = 10, imgPrint = TRUE,
                                            dir = paste0(result_dir))
    pathway_result_KEGG <- pathwayAnalysis_CP(geneList = important_entrez_ids,
                                              org = species[1], database = "KEGG",
                                              title = paste0(PC, "_Pathway_Results"),
                                              displayNum = 10, imgPrint = TRUE,
                                              dir = paste0(result_dir))
    
    ### enrichment object
    eo <- enrichGO(gene = important_entrez_ids, OrgDb = 'org.Mm.eg.db', readable = T, ont = "BP", pvalueCutoff = 0.05)
    eox <- setReadable(eo, 'org.Mm.eg.db', 'ENTREZID')
  }
  write.xlsx2(pathway_result_GO, file = paste0(result_dir, PC, "_GO_Pathway_Results_with_", length(important_genes), "_important_genes_", important_thresh, ".xlsx"),
              row.names = FALSE, sheetName = paste0("GO_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(result_dir, PC, "_KEGG_Pathway_Results_with_", length(important_genes), "_important_genes_", important_thresh, ".xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_Results"))
  
  ### dot plot
  p <- dotplot(eo, showCategory=10) + theme_classic(base_size = 40)
  p$theme$legend.key.size <- unit(2, "lines")
  ggsave(paste0(result_dir, "/DotPlot_", PC, "_Pathways_with_", length(important_genes), "_important_genes_", important_thresh, ".png"), plot = p, width = 22, height = 14, dpi = 300)
  
  ### upset plot
  p <- upsetplot(eo, n = 10) + theme(axis.text.y = element_text(size = 30))
  ggsave(paste0(result_dir, "/UpsetPlot_", PC, "_Pathways_with_", length(important_genes), "_important_genes_", important_thresh, ".png"), plot = p, width = 22, height = 12, dpi = 300)
  
  #
  ### GSEA
  #
  
  ### db preparation
  # MSIGDB
  if(species[1] == "human") {
    m_df <- msigdbr(species = "Homo sapiens") 
  } else if(species[1] == "mouse") {
    m_df <- msigdbr(species = "Mus musculus")
  }
  m_list <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  ### signature preparation
  # signat <- pca_contb[,PC]
  signat <- nrow(pca_contb):1
  names(signat) <- rownames(pca_contb)
  
  ### run GSEA
  GSEA_result <- run_gsea(gene_list = m_list, signature = list(signat), printPlot = FALSE)
  GSEA_result <- GSEA_result[order(GSEA_result$pval),]
  
  ### only get pathways that have pval < 0.001 & size > 30 & up-regulating results (enriched with important) only
  pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-03),
                                                      which(GSEA_result$size > 30)),
                                            which(abs(GSEA_result$NES) > 2))]
  if(length(pathways) < 5) {
    pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-03),
                                                        which(GSEA_result$size > 30)),
                                              which(abs(GSEA_result$NES) > 1.8))]
  }
  if(length(pathways) < 5) {
    pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-03),
                                                        which(GSEA_result$size > 30)),
                                              which(abs(GSEA_result$NES) > 1.6))]
  }
  
  ### run GSEA again with the significant result - plot printing
  result_dir2 <- paste0(result_dir, "GSEA/")
  dir.create(result_dir2, showWarnings = FALSE, recursive = TRUE)
  GSEA_result2 <- run_gsea(gene_list = m_list[pathways], signature = list(signat),
                           printPlot = TRUE, printPath = result_dir2)
  GSEA_result2 <- GSEA_result2[order(GSEA_result2$padj, GSEA_result2$pval),]
  
  ### write out the result file
  write.xlsx2(GSEA_result2, file = paste0(result_dir2, PC, "_Genes_GSEA_Results_msigdb.xlsx"),
              sheetName = "GSEA_Result", row.names = FALSE)
  
  ### same analyses with the given comparison
  if(!is.null(PC_Val)) {
    ### check rownames in pca map and meta.data and colnames in raw counts are the same
    if(identical(rownames(pca_map),
                 rownames(local_Seurat_Obj@meta.data)) && identical(rownames(local_Seurat_Obj@meta.data),
                                                                    colnames(local_Seurat_Obj@assays$RNA@counts))) {
      
      ### set two groups
      grp1_idx <- which(pca_map[,PC] >= PC_Val)
      grp2_idx <- which(pca_map[,PC] < PC_Val)
      
      ### set group info to the meta data
      local_Seurat_Obj@meta.data$PC_Group <- NA
      local_Seurat_Obj@meta.data$PC_Group[grp1_idx] <- "grp1"
      local_Seurat_Obj@meta.data$PC_Group[grp2_idx] <- "grp2"
      
      ### set the ident of the object with the group info
      local_Seurat_Obj <- SetIdent(object = local_Seurat_Obj,
                                   cells = rownames(local_Seurat_Obj@meta.data),
                                   value = local_Seurat_Obj@meta.data$PC_Group)
      
      ### get DE genes between two groups
      de_result <- FindMarkers(object = local_Seurat_Obj, test.use = "DESeq2",
                               ident.1 = "grp1", ident.2 = "grp2",
                               logfc.threshold = 0, min.pct = 0.1)
      
      ### order the DE reuslt
      de_result <- de_result[order(de_result$p_val_adj),]
      
      ### add pct for each library
      unique_devel <- unique(local_Seurat_Obj@meta.data[,target_col])
      for(devel in unique_devel) {
        ### grp1 - devel indicies
        devel_idx1 <- intersect(grp1_idx, which(local_Seurat_Obj@meta.data$Development == devel))
        
        ### add the columns
        de_result <- cbind(de_result, sapply(rownames(de_result), function(x) {
          r <- length(which(local_Seurat_Obj@assays$RNA@counts[x,devel_idx1] > 0)) / length(grp1_idx)
          return(round(r, digits = 3))
        }))
        
        ### change column name
        colnames(de_result)[ncol(de_result)] <- paste0("pct.1_", devel)
      }
      for(devel in unique_devel) {
        ### grp2 - devel indicies
        devel_idx2 <- intersect(grp2_idx, which(local_Seurat_Obj@meta.data$Development == devel))
        
        ### add the columns
        de_result <- cbind(de_result, sapply(rownames(de_result), function(x) {
          r <- length(which(local_Seurat_Obj@assays$RNA@counts[x,devel_idx2] > 0)) / length(grp2_idx)
          return(round(r, digits = 3))
        }))
        
        ### change column name
        colnames(de_result)[ncol(de_result)] <- paste0("pct.2_", devel)
      }
      
      ### add logFC for each library
      for(devel in unique_devel) {
        ### grp1 - devel indicies
        devel_idx1 <- intersect(grp1_idx, which(local_Seurat_Obj@meta.data$Development == devel))
        
        ### grp2 - devel indicies
        devel_idx2 <- intersect(grp2_idx, which(local_Seurat_Obj@meta.data$Development == devel))
        
        ### only if there are cells from the both group
        if(length(devel_idx1) > 0 && length(devel_idx2) > 0) {
          de_result <- cbind(de_result, sapply(rownames(de_result), function(x) {
            r <- log2(mean(local_Seurat_Obj@assays$RNA@data[x,devel_idx1]) / mean(local_Seurat_Obj@assays$RNA@data[x,devel_idx2]))
            return(round(r, digits = 8))
          }))  
          
          ### change column name
          colnames(de_result)[ncol(de_result)] <- paste0("logFC_", devel)
        }
      }
      
      ### new directory for DE
      result_dir2 <- paste0(result_dir, "PC_Comparison/", PC, "_", PC_Val, "/")
      dir.create(result_dir2, showWarnings = FALSE, recursive = TRUE)
      
      ### write out the DE result
      write.xlsx2(data.frame(Gene_Symbol=rownames(de_result),
                             de_result,
                             stringsAsFactors = FALSE, check.names = FALSE),
                  file = paste0(result_dir2, "DE_Result_", PC, "_", PC_Val, ".xlsx"),
                  sheetName = "DE_Result",
                  row.names = FALSE)
      
      ### at least there are 10 DE genes for the further analyses
      if(length(which(de_result$p_val_adj < 0.01)) > 10) {
        ### get important DE genes
        de_genes <- rownames(de_result)[which(de_result$p_val_adj < 0.01)]
        up_de_genes <- rownames(de_result)[intersect(which(de_result$p_val_adj < 0.01),
                                                     which(de_result$avg_log2FC > 0))]
        down_de_genes <- rownames(de_result)[intersect(which(de_result$p_val_adj < 0.01),
                                                       which(de_result$avg_log2FC < 0))]
        
        ### pathway analysis
        if(species[1] == "human") {
          if(one_side_only) {
            ### up genes
            if(length(up_de_genes) > 5) {
              ### get entrez ids for the genes
              de_entrez_ids <- mapIds(org.Hs.eg.db,
                                      up_de_genes,
                                      "ENTREZID", "SYMBOL")
              de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
              
              ### GO & KEGG
              pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                                      org = species[1], database = "GO",
                                                      title = paste0(PC, "_", PC_Val, "_Pathways_with_Up_Genes"),
                                                      displayNum = 10, imgPrint = TRUE,
                                                      dir = paste0(result_dir2))
              pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                                        org = species[1], database = "KEGG",
                                                        title = paste0(PC, "_", PC_Val, "_Pathways_with_Up_Genes"),
                                                        displayNum = 10, imgPrint = TRUE,
                                                        dir = paste0(result_dir2))
              write.xlsx2(pathway_result_GO, file = paste0(result_dir2, PC, "_", PC_Val, "_GO_Pathways_with_Up_", length(up_de_genes), "_DE_genes.xlsx"),
                          row.names = FALSE, sheetName = paste0("GO_Results"))
              write.xlsx2(pathway_result_KEGG, file = paste0(result_dir2, PC, "_", PC_Val, "_KEGG_Pathways_with_Up_", length(up_de_genes), "_DE_genes.xlsx"),
                          row.names = FALSE, sheetName = paste0("KEGG_Results"))
              
              ### enrichment object
              eo1 <- enrichGO(gene = de_entrez_ids, OrgDb = 'org.Hs.eg.db', readable = T, ont = "BP", pvalueCutoff = 0.05)
              eox1 <- setReadable(eo1, 'org.Hs.eg.db', 'ENTREZID')
              
              ### dot plot
              p <- dotplot(eo1, showCategory=10) + theme_classic(base_size = 40)
              p$theme$legend.key.size <- unit(2, "lines")
              ggsave(paste0(result_dir2, "/DotPlot_", PC, "_", PC_Val, "_GO_Pathway_Results_with_Up_", length(up_de_genes), "_DE_genes.png"), plot = p, width = 22, height = 14, dpi = 300)
              
              ### upset plot
              p <- upsetplot(eo1, n = 10) + theme(axis.text.y = element_text(size = 30))
              ggsave(paste0(result_dir2, "/UpsetPlot_", PC, "_", PC_Val, "_GO_Pathway_Results_with_Up_", length(up_de_genes), "_DE_genes.png"), plot = p, width = 20, height = 12, dpi = 300)
            }
            
            ### down genes
            if(length(down_de_genes) > 5) {
              ### get entrez ids for the genes
              de_entrez_ids <- mapIds(org.Hs.eg.db,
                                      down_de_genes,
                                      "ENTREZID", "SYMBOL")
              de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
              
              ### GO & KEGG
              pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                                      org = species[1], database = "GO",
                                                      title = paste0(PC, "_", PC_Val, "_Pathways_with_Down_Genes"),
                                                      displayNum = 10, imgPrint = TRUE,
                                                      dir = paste0(result_dir2))
              pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                                        org = species[1], database = "KEGG",
                                                        title = paste0(PC, "_", PC_Val, "_Pathways_with_Down_Genes"),
                                                        displayNum = 10, imgPrint = TRUE,
                                                        dir = paste0(result_dir2))
              write.xlsx2(pathway_result_GO, file = paste0(result_dir2, PC, "_", PC_Val, "_GO_Pathway_Results_with_Down_", length(down_de_genes), "_DE_genes.xlsx"),
                          row.names = FALSE, sheetName = paste0("GO_Results"))
              write.xlsx2(pathway_result_KEGG, file = paste0(result_dir2, PC, "_", PC_Val, "_KEGG_Pathway_Results_with_Down_", length(down_de_genes), "_DE_genes.xlsx"),
                          row.names = FALSE, sheetName = paste0("KEGG_Results"))
              
              ### enrichment object
              eo2 <- enrichGO(gene = de_entrez_ids, OrgDb = 'org.Hs.eg.db', readable = T, ont = "BP", pvalueCutoff = 0.05)
              eox2 <- setReadable(eo2, 'org.Hs.eg.db', 'ENTREZID')
              
              ### dot plot
              p <- dotplot(eo2, showCategory=10) + theme_classic(base_size = 40)
              p$theme$legend.key.size <- unit(2, "lines")
              ggsave(paste0(result_dir2, "/DotPlot_", PC, "_", PC_Val, "_GO_Pathway_Results_with_Down_", length(down_de_genes), "_DE_genes.png"), plot = p, width = 22, height = 14, dpi = 300)
              
              ### upset plot
              p <- upsetplot(eo2, n = 10) + theme(axis.text.y = element_text(size = 30))
              ggsave(paste0(result_dir2, "/UpsetPlot_", PC, "_", PC_Val, "_GO_Pathway_Results_with_Down_", length(down_de_genes), "_DE_genes.png"), plot = p, width = 22, height = 12, dpi = 300)
            }
          } else {
            ### get entrez ids for the genes
            de_entrez_ids <- mapIds(org.Hs.eg.db,
                                    de_genes,
                                    "ENTREZID", "SYMBOL")
            de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
            
            ### GO & KEGG
            pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                                    org = species[1], database = "GO",
                                                    title = paste0(PC, "_", PC_Val, "_Pathways_with_DE_Genes"),
                                                    displayNum = 10, imgPrint = TRUE,
                                                    dir = paste0(result_dir2))
            pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                                      org = species[1], database = "KEGG",
                                                      title = paste0(PC, "_", PC_Val, "_Pathways_DE_Genes"),
                                                      displayNum = 10, imgPrint = TRUE,
                                                      dir = paste0(result_dir2))
            write.xlsx2(pathway_result_GO, file = paste0(result_dir2, PC, "_", PC_Val, "_GO_Pathway_Results_with_", length(de_genes), "_DE_genes.xlsx"),
                        row.names = FALSE, sheetName = paste0("GO_Results"))
            write.xlsx2(pathway_result_KEGG, file = paste0(result_dir2, PC, "_", PC_Val, "_KEGG_Pathway_Results_with_", length(de_genes), "_DE_genes.xlsx"),
                        row.names = FALSE, sheetName = paste0("KEGG_Results"))
            
            ### enrichment object
            eo <- enrichGO(gene = de_entrez_ids, OrgDb = 'org.Hs.eg.db', readable = T, ont = "BP", pvalueCutoff = 0.05)
            eox <- setReadable(eo, 'org.Hs.eg.db', 'ENTREZID')
            
            ### dot plot
            p <- dotplot(eo, showCategory=10) + theme_classic(base_size = 40)
            p$theme$legend.key.size <- unit(2, "lines")
            ggsave(paste0(result_dir2, "/DotPlot_", PC, "_", PC_Val, "_GO_Pathway_Results_with_", length(de_genes), "_DE_genes.png"), plot = p, width = 22, height = 14, dpi = 300)
            
            ### upset plot
            p <- upsetplot(eo, n = 10) + theme(axis.text.y = element_text(size = 30))
            ggsave(paste0(result_dir2, "/UpsetPlot_", PC, "_", PC_Val, "_GO_Pathway_Results_with_", length(de_genes), "_DE_genes.png"), plot = p, width = 22, height = 12, dpi = 300)
          }
        } else if(species[1] == "mouse") {
          if(one_side_only) {
            ### up genes
            if(length(up_de_genes) > 5) {
              ### get entrez ids for the genes
              de_entrez_ids <- mapIds(org.Mm.eg.db,
                                      up_de_genes,
                                      "ENTREZID", "SYMBOL")
              de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
              
              ### GO & KEGG
              pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                                      org = species[1], database = "GO",
                                                      title = paste0(PC, "_", PC_Val, "_Pathways_with_Up_Genes"),
                                                      displayNum = 10, imgPrint = TRUE,
                                                      dir = paste0(result_dir2))
              pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                                        org = species[1], database = "KEGG",
                                                        title = paste0(PC, "_", PC_Val, "_Pathways_with_Up_Genes"),
                                                        displayNum = 10, imgPrint = TRUE,
                                                        dir = paste0(result_dir2))
              write.xlsx2(pathway_result_GO, file = paste0(result_dir2, PC, "_", PC_Val, "_GO_Pathway_Results_with_Up_", length(up_de_genes), "_DE_genes.xlsx"),
                          row.names = FALSE, sheetName = paste0("GO_Results"))
              write.xlsx2(pathway_result_KEGG, file = paste0(result_dir2, PC, "_", PC_Val, "_KEGG_Pathway_Results_with_Up_", length(up_de_genes), "_DE_genes.xlsx"),
                          row.names = FALSE, sheetName = paste0("KEGG_Results"))
              
              ### enrichment object
              eo1 <- enrichGO(gene = de_entrez_ids, OrgDb = 'org.Mm.eg.db', readable = T, ont = "BP", pvalueCutoff = 0.05)
              eox1 <- setReadable(eo1, 'org.Mm.eg.db', 'ENTREZID')
              
              ### dot plot
              p <- dotplot(eo1, showCategory=10) + theme_classic(base_size = 40)
              p$theme$legend.key.size <- unit(2, "lines")
              ggsave(paste0(result_dir2, "/DotPlot_", PC, "_", PC_Val, "_GO_Pathway_Results_with_Up_", length(up_de_genes), "_DE_genes.png"), plot = p, width = 22, height = 14, dpi = 300)
              
              ### upset plot
              p <- upsetplot(eo1, n = 10) + theme(axis.text.y = element_text(size = 30))
              ggsave(paste0(result_dir2, "/UpsetPlot_", PC, "_", PC_Val, "_GO_Pathway_Results_with_Up_", length(up_de_genes), "_DE_genes.png"), plot = p, width = 22, height = 12, dpi = 300)
            }
            
            ### down genes
            if(length(down_de_genes) > 5) {
              ### get entrez ids for the genes
              de_entrez_ids <- mapIds(org.Mm.eg.db,
                                      down_de_genes,
                                      "ENTREZID", "SYMBOL")
              de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
              
              ### GO & KEGG
              pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                                      org = species[1], database = "GO",
                                                      title = paste0(PC, "_", PC_Val, "_Pathways_with_Down_Genes"),
                                                      displayNum = 10, imgPrint = TRUE,
                                                      dir = paste0(result_dir2))
              pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                                        org = species[1], database = "KEGG",
                                                        title = paste0(PC, "_", PC_Val, "_Pathways_with_Down_Genes"),
                                                        displayNum = 10, imgPrint = TRUE,
                                                        dir = paste0(result_dir2))
              write.xlsx2(pathway_result_GO, file = paste0(result_dir2, PC, "_", PC_Val, "_GO_Pathway_Results_with_Down_", length(down_de_genes), "_DE_genes.xlsx"),
                          row.names = FALSE, sheetName = paste0("GO_Results"))
              write.xlsx2(pathway_result_KEGG, file = paste0(result_dir2, PC, "_", PC_Val, "_KEGG_Pathway_Results_with_Down_", length(down_de_genes), "_DE_genes.xlsx"),
                          row.names = FALSE, sheetName = paste0("KEGG_Results"))
              
              ### enrichment object
              eo2 <- enrichGO(gene = de_entrez_ids, OrgDb = 'org.Mm.eg.db', readable = T, ont = "BP", pvalueCutoff = 0.05)
              eox2 <- setReadable(eo2, 'org.Mm.eg.db', 'ENTREZID')
              
              ### dot plot
              p <- dotplot(eo2, showCategory=10) + theme_classic(base_size = 40)
              p$theme$legend.key.size <- unit(2, "lines")
              ggsave(paste0(result_dir2, "/DotPlot_", PC, "_", PC_Val, "_GO_Pathway_Results_with_Down_", length(down_de_genes), "_DE_genes.png"), plot = p, width = 22, height = 14, dpi = 300)
              
              ### upset plot
              p <- upsetplot(eo2, n = 10) + theme(axis.text.y = element_text(size = 30))
              ggsave(paste0(result_dir2, "/UpsetPlot_", PC, "_", PC_Val, "_GO_Pathway_Results_with_Down_", length(down_de_genes), "_DE_genes.png"), plot = p, width = 22, height = 12, dpi = 300)
            }
          } else {
            ### get entrez ids for the genes
            de_entrez_ids <- mapIds(org.Mm.eg.db,
                                    de_genes,
                                    "ENTREZID", "SYMBOL")
            de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
            
            ### GO & KEGG
            pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                                    org = species[1], database = "GO",
                                                    title = paste0(PC, "_", PC_Val, "_Pathways_with_DE_Genes"),
                                                    displayNum = 10, imgPrint = TRUE,
                                                    dir = paste0(result_dir2))
            pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                                      org = species[1], database = "KEGG",
                                                      title = paste0(PC, "_", PC_Val, "_Pathways_with_DE_Genes"),
                                                      displayNum = 10, imgPrint = TRUE,
                                                      dir = paste0(result_dir2))
            write.xlsx2(pathway_result_GO, file = paste0(result_dir2, PC, "_", PC_Val, "_GO_Pathway_Results_with_", length(de_genes), "_DE_genes.xlsx"),
                        row.names = FALSE, sheetName = paste0("GO_Results"))
            write.xlsx2(pathway_result_KEGG, file = paste0(result_dir2, PC, "_", PC_Val, "_KEGG_Pathway_Results_with_", length(de_genes), "_DE_genes.xlsx"),
                        row.names = FALSE, sheetName = paste0("KEGG_Results"))
            
            ### enrichment object
            eo <- enrichGO(gene = de_entrez_ids, OrgDb = 'org.Mm.eg.db', readable = T, ont = "BP", pvalueCutoff = 0.05)
            eox <- setReadable(eo, 'org.Mm.eg.db', 'ENTREZID')
            
            ### dot plot
            p <- dotplot(eo, showCategory=10) + theme_classic(base_size = 40)
            p$theme$legend.key.size <- unit(2, "lines")
            ggsave(paste0(result_dir2, "/DotPlot_", PC, "_", PC_Val, "_GO_Pathway_Results_with_", length(de_genes), "_DE_genes.png"), plot = p, width = 22, height = 14, dpi = 300)
            
            ### upset plot
            p <- upsetplot(eo, n = 10) + theme(axis.text.y = element_text(size = 30))
            ggsave(paste0(result_dir2, "/UpsetPlot_", PC, "_", PC_Val, "_GO_Pathway_Results_with_", length(de_genes), "_DE_genes.png"), plot = p, width = 22, height = 12, dpi = 300)
          }
        }
        
        ### GSEA
        ### signature preparation
        signat <- de_result$avg_log2FC
        names(signat) <- rownames(de_result)
        
        ### run GSEA
        GSEA_result <- run_gsea(gene_list = m_list, signature = list(signat), printPlot = FALSE)
        GSEA_result <- GSEA_result[order(GSEA_result$pval),]
        
        ### only get pathways that have pval < 0.001 & size > 30 & up-regulating results (enriched with important) only
        pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-03),
                                                            which(GSEA_result$size > 30)),
                                                  which(abs(GSEA_result$NES) > 2))]
        
        if(length(pathways) < 5) {
          pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-03),
                                                              which(GSEA_result$size > 30)),
                                                    which(abs(GSEA_result$NES) > 1.8))]
        }
        if(length(pathways) < 5) {
          pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-03),
                                                              which(GSEA_result$size > 30)),
                                                    which(abs(GSEA_result$NES) > 1.6))]
        }
        
        ### run GSEA again with the significant result - plot printing
        result_dir2 <- paste0(result_dir2, "GSEA/")
        dir.create(result_dir2, showWarnings = FALSE, recursive = TRUE)
        GSEA_result2 <- run_gsea(gene_list = m_list[pathways], signature = list(signat),
                                 printPlot = TRUE, printPath = result_dir2)
        GSEA_result2 <- GSEA_result2[order(GSEA_result2$padj, GSEA_result2$pval),]
        
        ### write out the result file
        write.xlsx2(GSEA_result2, file = paste0(result_dir2, PC, "_", PC_Val, "_DE_Genes_GSEA_Results_msigdb.xlsx"),
                    sheetName = "GSEA_Result", row.names = FALSE)
      }
      
    } else {
      stop("ERROR: multiple_analyses_in_one() - unidentical row or col names.")
    }
  }
  
}

### load the Seurat object and save the object name
tmp_env <- new.env()
load(paste0(Robj1_path), tmp_env)
obj_name <- ls(tmp_env)
assign("Combined_Seurat_Obj", get(obj_name, envir = tmp_env))
rm(tmp_env)
gc()

### run PCA
Combined_Seurat_Obj <- FindVariableFeatures(Combined_Seurat_Obj)
Combined_Seurat_Obj <- ScaleData(Combined_Seurat_Obj)
Combined_Seurat_Obj <- RunPCA(Combined_Seurat_Obj, npcs = 10)

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

### set the ident of the object with the HSPC type
Combined_Seurat_Obj <- SetIdent(object = Combined_Seurat_Obj,
                                cells = rownames(Combined_Seurat_Obj@meta.data),
                                value = Combined_Seurat_Obj@meta.data$Development)

### check whether the orders are the same
print(identical(names(Idents(object = Combined_Seurat_Obj)), rownames(Combined_Seurat_Obj@meta.data)))

### rownames in the meta.data should be in the same order as colnames in the counts
Combined_Seurat_Obj@meta.data <- Combined_Seurat_Obj@meta.data[colnames(Combined_Seurat_Obj@assays$RNA@counts),]

### active assay = "RNA"
Combined_Seurat_Obj@active.assay <- "RNA"

### get unique idents
unique_idents <- unique(Combined_Seurat_Obj@meta.data$Development)

### a function for color brewer
cell_pal <- function(cell_vars, pal_fun) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories)), categories)
    return(pal[cell_vars])
  }
}


#' @title Plot Slingshot output
#' @name plot-SlingshotDataSet
#' @aliases plot-SlingshotDataSet plot,SlingshotDataSet,ANY-method
#'
#' @description Tools for visualizing lineages inferred by \code{slingshot}.
#'
#' @param x a \code{SlingshotDataSet} with results to be plotted.
#' @param type character, the type of output to be plotted, can be one of
#'   \code{"lineages"}, \code{"curves"}, or \code{"both"} (by partial matching),
#'   see Details for more.
#' @param linInd integer, an index indicating which lineages should be plotted
#'   (default is to plot all lineages). If \code{col} is a vector, it will be
#'   subsetted by \code{linInd}.
#' @param show.constraints logical, whether or not the user-specified initial
#'   and terminal clusters should be specially denoted by green and red dots,
#'   respectively.
#' @param add logical, indicates whether the output should be added to an
#'   existing plot.
#' @param dims numeric, which dimensions to plot (default is \code{1:2}).
#' @param asp numeric, the y/x aspect ratio, see \code{\link{plot.window}}.
#' @param cex numeric, amount by which points should be magnified, see
#'   \code{\link{par}}.
#' @param lwd numeric, the line width, see \code{\link{par}}.
#' @param col character or numeric, color(s) for lines, see \code{\link{par}}.
#' @param ... additional parameters to be passed to \code{\link{lines}}.
#'
#' @details If \code{type == 'lineages'}, straight line connectors between
#'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
#'   principal curves will be plotted.
#'
#' @details When \code{type} is not specified, the function will first check the
#'   \code{curves} slot and plot the curves, if present. Otherwise,
#'   \code{lineages} will be plotted, if present.
#'
#' @return returns \code{NULL}.
#'
#' @examples
#' data("slingshotExample")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' sds <- slingshot(rd, cl, start.clus = "1")
#' plot(sds, type = 'b')
#'
#' # add to existing plot
#' plot(rd, col = 'grey50')
#' lines(sds, lwd = 3)
#'
#' @import graphics
#' @import grDevices
#' @export
setMethod(
  f = "plot",
  signature = signature(x = "SlingshotDataSet"),
  definition = function(x, type = NULL,
                        linInd = NULL,
                        show.constraints = FALSE,
                        constraints.col = NULL,
                        add = FALSE,
                        dims = seq_len(2),
                        asp = 1,
                        cex = 2,
                        lwd = 2,
                        col = 1,
                        ...) {
    col <- rep(col, length(slingLineages(x)))
    curves <- FALSE
    lineages <- FALSE
    if(is.null(type)){
      if(length(slingCurves(x)) > 0){
        type <- 'curves'
      }else if(length(slingLineages(x)) > 0){
        type <- 'lineages'
      }else{
        stop('No lineages or curves detected.')
      }
    }else{
      type <- c('curves','lineages','both')[pmatch(type,
                                                   c('curves','lineages','both'))]
      if(is.na(type)){
        stop('Unrecognized type argument.')
      }
    }
    
    if(type %in% c('lineages','both')){
      lineages <- TRUE
    }
    if(type %in% c('curves','both')){
      curves <- TRUE
    }
    
    if(lineages & (length(slingLineages(x))==0)){
      stop('No lineages detected.')
    }
    if(curves & (length(slingCurves(x))==0)){
      stop('No curves detected.')
    }
    
    if(is.null(linInd)){
      linInd <- seq_along(slingLineages(x))
    }else{
      linInd <- as.integer(linInd)
      if(!all(linInd %in% seq_along(slingLineages(x)))){
        if(any(linInd %in% seq_along(slingLineages(x)))){
          linInd.removed <-
            linInd[! linInd %in% seq_along(slingLineages(x))]
          linInd <-
            linInd[linInd %in% seq_along(slingLineages(x))]
          message('Unrecognized lineage indices (linInd): ',
                  paste(linInd.removed, collapse = ", "))
        }else{
          stop('None of the provided lineage indices',
               ' (linInd) were found.')
        }
      }
    }
    
    if(lineages){
      X <- reducedDim(x)
      clusterLabels <- slingClusterLabels(x)
      connectivity <- slingAdjacency(x)
      clusters <- rownames(connectivity)
      nclus <- nrow(connectivity)
      centers <- t(vapply(clusters,function(clID){
        w <- clusterLabels[,clID]
        return(apply(X, 2, weighted.mean, w = w))
      }, rep(0,ncol(X))))
      rownames(centers) <- clusters
      X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
      clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, ,
                                     drop = FALSE]
      linC <- slingParams(x)
      clus2include <- unique(unlist(slingLineages(x)[linInd]))
    }
    
    if(!add){
      xs <- NULL
      ys <- NULL
      if(lineages){
        xs <- c(xs, centers[,dims[1]])
        ys <- c(ys, centers[,dims[2]])
      }
      if(curves){
        npoints <- nrow(slingCurves(x)[[1]]$s)
        xs <- c(xs, as.numeric(vapply(slingCurves(x),
                                      function(c){ c$s[,dims[1]] }, rep(0,npoints))))
        ys <- c(ys, as.numeric(vapply(slingCurves(x),
                                      function(c){ c$s[,dims[2]] }, rep(0,npoints))))
      }
      plot(x = NULL, y = NULL, asp = asp,
           xlim = range(xs), ylim = range(ys),
           xlab = colnames(reducedDim(x))[dims[1]],
           ylab = colnames(reducedDim(x))[dims[2]])
    }
    
    if(lineages){
      for(i in seq_len(nclus-1)){
        for(j in seq(i+1,nclus)){
          if(connectivity[i,j]==1 &
             all(clusters[c(i,j)] %in% clus2include)){
            lines(centers[c(i,j), dims],
                  lwd = lwd, col = col[1], ...)
          }
        }
      }
      points(centers[clusters %in% clus2include, dims],
             cex = cex+1, pch = 16, col = col[1])
      if(show.constraints && !is.null(constraints.col)){
        for(const in names(constraints.col)) {
          points(centers[clusters %in% const, dims,
                         drop=FALSE], cex = cex-0.5,
                 col = constraints.col[const], pch = 16)
          # text(x = centers[clusters %in% const, dims[1]]+0,
          #      y = centers[clusters %in% const, dims[2]]+2,
          #      labels = const,
          #      font = 2,
          #      cex = cex-0.5,
          #      col = "black")
        }
      }
    }
    if(curves){
      for(ii in seq_along(slingCurves(x))[linInd]){
        c <- slingCurves(x)[[ii]]
        lines(c$s[c$ord, dims], lwd = lwd, col = col[ii], ...)
      }
    }
    invisible(NULL)
  }
)

### remove ADULT cells
subset_Seurat_Obj <- subset(Combined_Seurat_Obj, idents = unique_idents[-which(unique_idents == "ADULT")])

### rownames in the meta.data should be in the same order as colnames in the counts
subset_Seurat_Obj@assays$RNA@counts <- subset_Seurat_Obj@assays$RNA@counts[,rownames(subset_Seurat_Obj@meta.data)]

### run PCA
subset_Seurat_Obj <- FindVariableFeatures(subset_Seurat_Obj)
subset_Seurat_Obj <- ScaleData(subset_Seurat_Obj)
subset_Seurat_Obj <- RunPCA(subset_Seurat_Obj, npcs = 15)


### set the ident of the object with the HSPC type
subset_Seurat_Obj <- SetIdent(object = subset_Seurat_Obj,
                              cells = rownames(subset_Seurat_Obj@meta.data),
                              value = subset_Seurat_Obj@meta.data$HSPC)

### check whether the orders are the same
print(identical(names(Idents(object = subset_Seurat_Obj)), rownames(subset_Seurat_Obj@meta.data)))

### new output directory
type <- "MPP2"
outputDir2 <- paste0(outputDir, type, "/")
dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)

### split the Seurat obj based on HSPC info
subset_Seurat_Obj2 <- subset(subset_Seurat_Obj, idents=type)

### order the meta data by developmental time
subset_Seurat_Obj2@meta.data <- subset_Seurat_Obj2@meta.data[order(factor(subset_Seurat_Obj2@meta.data$Development,
                                                                          levels = time_points)),]

### rownames in the meta.data should be in the same order as colnames in the counts
subset_Seurat_Obj2@assays$RNA@counts <- subset_Seurat_Obj2@assays$RNA@counts[,rownames(subset_Seurat_Obj2@meta.data)]

### run PCA
subset_Seurat_Obj2 <- RunPCA(subset_Seurat_Obj2, npcs = 10)
pca_map <- Embeddings(subset_Seurat_Obj2, reduction = "pca")[rownames(subset_Seurat_Obj2@meta.data),1:10]

### get slingshot object
slingshot_obj <- slingshot(pca_map,
                           clusterLabels = subset_Seurat_Obj2@meta.data$Development, 
                           reducedDim = "PCA")

### get colors for the clustering result
cell_colors_clust <- cell_pal(unique(subset_Seurat_Obj2@meta.data$Development), hue_pal())

### Trajectory inference
png(paste0(outputDir2, "Trajectory_Inference_Without_Adult_PCA.png"), width = 1800, height = 1200, res = 300)
plot(reducedDim(slingshot_obj),
     main=paste(type, "Trajectory Inference Without Adult"),
     col = cell_colors_clust[subset_Seurat_Obj2@meta.data$Development],
     pch = 19, cex = 1)
lines(slingshot_obj, lwd = 3, type = "lineages", col = "black",
      show.constraints = TRUE, constraints.col = cell_colors_clust)
legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
       pch = 19)
dev.off()

### Construct a monocle cds
monocle_cds <- newCellDataSet(as(as.matrix(subset_Seurat_Obj2@assays$RNA@data), 'sparseMatrix'),
                              phenoData = new('AnnotatedDataFrame', data = subset_Seurat_Obj2@meta.data),
                              featureData = new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(subset_Seurat_Obj2@assays$RNA@data),
                                                                                        row.names = row.names(subset_Seurat_Obj2@assays$RNA@data),
                                                                                        stringsAsFactors = FALSE, check.names = FALSE)),
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

### run monocle
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- reduceDimension(monocle_cds, reduction_method = "DDRTree")
monocle_cds <- orderCells(monocle_cds)

### draw monocle plots
p <- plot_cell_trajectory(monocle_cds, color_by = "Development", cell_size = 5, cell_link_size = 3, show_branch_points = FALSE) +
  labs(color="") +
  theme_classic(base_size = 36) +
  theme(legend.position = "top",
        legend.title = element_text(size = 36),
        legend.text = element_text(size = 30))
ggsave(file = paste0(outputDir2, "Trajectory_Inference_Without_Adult_Monocle2.png"),
       plot = p,
       width = 15, height = 10, dpi = 350)

p <- plot_complex_cell_trajectory(monocle_cds, color_by = "Development", cell_size = 5, cell_link_size = 3, show_branch_points = FALSE) +
  labs(color="") +
  theme_classic(base_size = 36) +
  theme(legend.position = "top",
        legend.title = element_text(size = 36),
        legend.text = element_text(size = 30))
ggsave(file = paste0(outputDir2, "Trajectory_Inference_Without_Adult_Complex_Monocle2.png"),
       plot = p,
       width = 15, height = 10, dpi = 350)

### DE & pathway analyses
#
# PC1
multiple_analyses_in_one(Seurat_Object = subset_Seurat_Obj2,
                         target_ident = NULL,
                         target_col = "Development",
                         target_col_factor_level = unique(subset_Seurat_Obj2@meta.data$Development),
                         one_side_only = TRUE,
                         species = "mouse",
                         PC = "PC_1",
                         PC_Val = -10,
                         important_thresh = 0.1,
                         garbage_thresh = 1e-04,
                         result_dir = paste0(outputDir2, "PC1_-10/"))



### new output directory
type <- "LTHSC"
outputDir2 <- paste0(outputDir, type, "/")
dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)

### split the Seurat obj based on HSPC info
subset_Seurat_Obj2 <- subset(subset_Seurat_Obj, idents=type)

### order the meta data by developmental time
subset_Seurat_Obj2@meta.data <- subset_Seurat_Obj2@meta.data[order(factor(subset_Seurat_Obj2@meta.data$Development,
                                                                          levels = time_points)),]

### rownames in the meta.data should be in the same order as colnames in the counts
subset_Seurat_Obj2@assays$RNA@counts <- subset_Seurat_Obj2@assays$RNA@counts[,rownames(subset_Seurat_Obj2@meta.data)]

### run PCA
subset_Seurat_Obj2 <- RunPCA(subset_Seurat_Obj2, npcs = 10)
pca_map <- Embeddings(subset_Seurat_Obj2, reduction = "pca")[rownames(subset_Seurat_Obj2@meta.data),1:10]

### get slingshot object
slingshot_obj <- slingshot(pca_map,
                           clusterLabels = subset_Seurat_Obj2@meta.data$Development, 
                           reducedDim = "PCA")

### get colors for the clustering result
cell_colors_clust <- cell_pal(unique(subset_Seurat_Obj2@meta.data$Development), hue_pal())

### Trajectory inference
png(paste0(outputDir2, "Trajectory_Inference_Without_Adult_PCA.png"), width = 1800, height = 1200, res = 300)
plot(reducedDim(slingshot_obj),
     main=paste(type, "Trajectory Inference Without Adult"),
     col = cell_colors_clust[subset_Seurat_Obj2@meta.data$Development],
     pch = 19, cex = 1)
lines(slingshot_obj, lwd = 3, type = "lineages", col = "black",
      show.constraints = TRUE, constraints.col = cell_colors_clust)
legend("topleft", legend = names(cell_colors_clust), col = cell_colors_clust,
       pch = 19)
dev.off()

### Construct a monocle cds
monocle_cds <- newCellDataSet(as(as.matrix(subset_Seurat_Obj2@assays$RNA@data), 'sparseMatrix'),
                              phenoData = new('AnnotatedDataFrame', data = subset_Seurat_Obj2@meta.data),
                              featureData = new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(subset_Seurat_Obj2@assays$RNA@data),
                                                                                        row.names = row.names(subset_Seurat_Obj2@assays$RNA@data),
                                                                                        stringsAsFactors = FALSE, check.names = FALSE)),
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

### run monocle
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- reduceDimension(monocle_cds, reduction_method = "DDRTree")
monocle_cds <- orderCells(monocle_cds)

### draw monocle plots
p <- plot_cell_trajectory(monocle_cds, color_by = "Development", cell_size = 5, cell_link_size = 3, show_branch_points = FALSE) +
  labs(color="") +
  theme_classic(base_size = 36) +
  theme(legend.position = "top",
        legend.title = element_text(size = 36),
        legend.text = element_text(size = 30))
ggsave(file = paste0(outputDir2, "Trajectory_Inference_Without_Adult_Monocle2.png"),
       plot = p,
       width = 15, height = 10, dpi = 350)

p <- plot_complex_cell_trajectory(monocle_cds, color_by = "Development", cell_size = 5, cell_link_size = 3, show_branch_points = FALSE) +
  labs(color="") +
  theme_classic(base_size = 36) +
  theme(legend.position = "top",
        legend.title = element_text(size = 36),
        legend.text = element_text(size = 30))
ggsave(file = paste0(outputDir2, "Trajectory_Inference_Without_Adult_Complex_Monocle2.png"),
       plot = p,
       width = 15, height = 10, dpi = 350)

### DE & pathway analyses
#
# PC1
multiple_analyses_in_one(Seurat_Object = subset_Seurat_Obj2,
                         target_ident = NULL,
                         target_col = "Development",
                         target_col_factor_level = unique(subset_Seurat_Obj2@meta.data$Development),
                         one_side_only = TRUE,
                         species = "mouse",
                         PC = "PC_1",
                         PC_Val = -30,
                         important_thresh = 0.1,
                         garbage_thresh = 1e-04,
                         result_dir = paste0(outputDir2, "PC1_-30/"))

### updated stroma RDS file
Updated_Seurat_Obj <- readRDS(file = Robj2_path)

### create new column for the analysis
Updated_Seurat_Obj@meta.data$Dev_Anno <- paste0(Updated_Seurat_Obj@meta.data$Development, "_",
                                                Updated_Seurat_Obj@meta.data$Annotation)

### determine necessary variables
time_points <- c("E16", "E18", "P0", "ADULT")
HSPC_populations <- c("LTHSC", "STHSC", "MPP2", "MPP3", "MPP4")

### set the ident of the object with Tissue type
Updated_Seurat_Obj <- SetIdent(object = Updated_Seurat_Obj,
                               cells = rownames(Updated_Seurat_Obj@meta.data),
                               value = Updated_Seurat_Obj@meta.data$Tissue)

### check whether the orders are the same
print(identical(names(Idents(object = Updated_Seurat_Obj)), rownames(Updated_Seurat_Obj@meta.data)))

### run UMAP
Updated_Seurat_Obj <- FindVariableFeatures(Updated_Seurat_Obj)
Updated_Seurat_Obj <- ScaleData(Updated_Seurat_Obj)
Updated_Seurat_Obj <- RunUMAP(Updated_Seurat_Obj, dims = 1:5)


#
### Additional on 08/23/2021
### Trent wants to find potential markers for flow
### DEG for each of the T cell clusters in the E16.5 Heme objects vs all other cells
#

### load libraries
if(!require(Seurat, quietly = TRUE)) {
  install.packages("Seurat")
  require(Seurat, quietly = TRUE)
}
if(!require(xlsx, quietly = TRUE)) {
  install.packages("xlsx")
  require(xlsx, quietly = TRUE)
}
if(!require(ggplot2, quietly = TRUE)) {
  install.packages("ggplot2")
  require(ggplot2, quietly = TRUE)
}

### load the object
load("./data/New_Objects/HemeE16_Final.Robj")

### active assay = "RNA"
### we are interested in DE genes of scRNA-Seq here, not CITE-Seq
HemeE16_regress@active.assay <- "RNA"

### Each TCell vs other non-Tcell groups
HemeE16_regress$Group_Anno <- as.character(HemeE16_regress$Final_Annotation)
HemeE16_regress$Group_Anno[!grepl(pattern = "TCell", HemeE16_regress$Group_Anno, fixed = TRUE)] <- "Others"

### set the ident of the object with the group info
HemeE16_regress <- SetIdent(object = HemeE16_regress,
                            cells = rownames(HemeE16_regress@meta.data),
                            value = HemeE16_regress$Group_Anno)

### set output directory
output_dir <- "./results/DE_For_Flow/"

### run DE analysis for each T cell cluster
de_result_list <- vector("list", length = 6)
for(i in 1:6) {
  
  ### get DE genes between two groups
  de_result <- FindMarkers(object = HemeE16_regress,
                           test.use = "wilcox",
                           ident.1 = paste0("TCell", i),
                           ident.2 = "Others",
                           logfc.threshold = 0.2,
                           min.pct = 0.2)
  
  ### order the DE reuslt
  de_result <- de_result[order(de_result$p_val_adj),]
  
  ### put the result in the list
  de_result_list[[i]] <- de_result
  
  ### data frame
  de_result <- data.frame(Gene_Symbol=rownames(de_result),
                          de_result,
                          stringsAsFactors = FALSE, check.names = FALSE)
  
  ### write out the DE result
  write.xlsx2(de_result,
              file = paste0(output_dir,
                            "TCell", i, "_vs_Others_DE_Result.xlsx"),
              sheetName = "DE_Result",
              row.names = FALSE)
  
  gc()
  
}

### were there any intersections?
### create a combined table
result_table <- de_result_list[[1]][which(de_result_list[[1]]$p_val_adj < 1e-20),]
result_table$Count <- 1
for(i in 2:6) {
  
  temp <- de_result_list[[i]][which(de_result_list[[i]]$p_val_adj < 1e-20),]
  temp$Count <- 1
  
  commons <- intersect(rownames(result_table), rownames(temp))
  new_ones <- setdiff(rownames(temp), rownames(result_table))
  
  if(length(commons) > 0) {
    for(j in 1:length(commons)) {
      for(k in 1:(ncol(result_table)-1)) {
        result_table[commons[j],k] <- paste(result_table[commons[j],k],
                                            temp[commons[j],k],
                                            sep = ";")
      }
      result_table[commons[j],"Count"] <- result_table[commons[j],"Count"] + temp[commons[j], "Count"]
    }
  }
  
  if(length(new_ones) > 0) {
    result_table <- rbind(result_table,
                          temp[new_ones,,drop=FALSE])
  }
  
}

### order it by appearance
result_table <- result_table[order(-result_table$Count),]

### save the result_table
write.xlsx2(data.frame(Genes=rownames(result_table),
                       result_table,
                       stringsAsFactors = FALSE, check.names = FALSE),
            file = paste0(output_dir,
                          "Combined_DE_Result_Table.xlsx"),
            sheetName = "DE_Result_Table",
            row.names = FALSE)

### How about all the Tcell vs the others?
HemeE16_regress$Group_Anno2 <- as.character(HemeE16_regress$Group_Anno)
HemeE16_regress$Group_Anno2[grep(pattern = "TCell", HemeE16_regress$Group_Anno2, fixed = TRUE)] <- "TCells"

### set the ident of the object with the group info
HemeE16_regress <- SetIdent(object = HemeE16_regress,
                            cells = rownames(HemeE16_regress@meta.data),
                            value = HemeE16_regress$Group_Anno2)

### get DE genes between two groups
de_result <- FindMarkers(object = HemeE16_regress,
                         test.use = "wilcox",
                         ident.1 = "TCells",
                         ident.2 = "Others",
                         logfc.threshold = 0.2,
                         min.pct = 0.2)

### order the DE reuslt
de_result <- de_result[order(de_result$p_val_adj),]

### data frame
de_result <- data.frame(Gene_Symbol=rownames(de_result),
                        de_result,
                        stringsAsFactors = FALSE, check.names = FALSE)

### write out the DE result
write.xlsx2(de_result,
            file = paste0(output_dir,
                          "All_TCell_vs_Others_DE_Result.xlsx"),
            sheetName = "DE_Result",
            row.names = FALSE)

### 09/09/2021
### 1. flow-friendly markers for the mast cell population at E16.5
### 2. additional markers that might distinguish between the large T cell cluster and the smaller T cell cluster (that is still present at P0)

### see how the clusters look like
DimPlot(HemeE16_regress, reduction = "umap", group.by = "Final_Annotation", pt.size = 2, label = TRUE)

### Mast cell cluster vs the other cells
HemeE16_regress$Group_Anno3 <- "Others"
HemeE16_regress$Group_Anno3[grep(pattern = "^Mast", HemeE16_regress$Final_Annotation, fixed = FALSE)] <- "Mast"

### see how the clusters look like again
DimPlot(HemeE16_regress, reduction = "umap", group.by = "Group_Anno3", pt.size = 2, label = TRUE)

### set the ident of the object with the group info
HemeE16_regress <- SetIdent(object = HemeE16_regress,
                            cells = rownames(HemeE16_regress@meta.data),
                            value = HemeE16_regress$Group_Anno3)

### get DE genes between two groups
de_result <- FindMarkers(object = HemeE16_regress,
                         test.use = "wilcox",
                         ident.1 = "Mast",
                         ident.2 = "Others",
                         logfc.threshold = 0.2,
                         min.pct = 0.2)

### order the DE reuslt
de_result <- de_result[order(de_result$p_val_adj),]

### data frame
de_result <- data.frame(Gene_Symbol=rownames(de_result),
                        de_result,
                        stringsAsFactors = FALSE, check.names = FALSE)

### write out the DE result
write.xlsx2(de_result,
            file = paste0(output_dir,
                          "All_Mast_vs_Others_DE_Result.xlsx"),
            sheetName = "DE_Result",
            row.names = FALSE)

### DE analysis between large t cell cluster vs small t cell cluster
### but we have to check the small clusters are the same between E16.5 and P0

### load P0 data
load("./data/New_Objects/HemeP0_Final.Robj")

### active assay = "RNA"
### we are interested in DE genes of scRNA-Seq here, not CITE-Seq
HemeP0_regress@active.assay <- "RNA"

### combine HemeE16 and HemeP0
HemeE16_regress$Final_Annotation2 <- paste0("Heme16_", HemeE16_regress$Final_Annotation)
HemeP0_regress$Final_Annotation2 <- paste0("P0_", HemeP0_regress$Final_Annotation)
Heme_E16_P0_Obj <- merge(x = HemeE16_regress, y = HemeP0_regress, add.cell.ids = c("E16", "P0"))

### normalization
Heme_E16_P0_Obj <- NormalizeData(Heme_E16_P0_Obj,
                                 normalization.method = "LogNormalize", scale.factor = 10000)

### find variable genes
Heme_E16_P0_Obj <- FindVariableFeatures(Heme_E16_P0_Obj,
                                        selection.method = "vst", nfeatures = 2000)

### scaling
Heme_E16_P0_Obj <- ScaleData(Heme_E16_P0_Obj,
                             vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))

### run pca & umap
Heme_E16_P0_Obj <- RunPCA(Heme_E16_P0_Obj,
                          features = VariableFeatures(object = Heme_E16_P0_Obj),
                          npcs = 15)
Heme_E16_P0_Obj <- RunUMAP(Heme_E16_P0_Obj, dims = 1:15)

### see whether TCell5 clusters of E16 and P0 are clustered together
p <- DimPlot(Heme_E16_P0_Obj, reduction = "umap", group.by = "Final_Annotation2", pt.size = 2, label = TRUE) +
  labs(title = "") +
  theme_classic(base_size = 30) +
  theme(legend.position = "none")
ggsave(file = paste0(output_dir, "UMAP_Heme_E16_P0.pdf"), plot = p,
       width = 15, height = 10, dpi = 350)

### E16 large t cell cluster vs small t cell cluster
HemeE16_regress$Group_Anno4 <- "Others"
HemeE16_regress$Group_Anno4[which(HemeE16_regress$Final_Annotation %in% c("TCell1", "TCell2", "TCell3", "TCell4", "TCell6"))] <- "Large_T"
HemeE16_regress$Group_Anno4[which(HemeE16_regress$Final_Annotation == "TCell5")] <- "Small_T"

### check the Group_Anno4
DimPlot(HemeE16_regress, reduction = "umap", group.by = "Group_Anno4", pt.size = 2, label = TRUE)

### set the ident of the object with the group info
HemeE16_regress <- SetIdent(object = HemeE16_regress,
                            cells = rownames(HemeE16_regress@meta.data),
                            value = HemeE16_regress$Group_Anno4)

### get DE genes between two groups
de_result <- FindMarkers(object = HemeE16_regress,
                         test.use = "wilcox",
                         ident.1 = "Large_T",
                         ident.2 = "Small_T",
                         logfc.threshold = 0.2,
                         min.pct = 0.2)

### order the DE reuslt
de_result <- de_result[order(de_result$p_val_adj),]

### data frame
de_result <- data.frame(Gene_Symbol=rownames(de_result),
                        de_result,
                        stringsAsFactors = FALSE, check.names = FALSE)

### write out the DE result
write.xlsx2(de_result,
            file = paste0(output_dir,
                          "Large_TCell_vs_Small_TCell_DE_Result.xlsx"),
            sheetName = "DE_Result",
            row.names = FALSE)


### How about all the Tcell vs the others?
HemeE16_regress$Group_Anno2 <- as.character(HemeE16_regress$Group_Anno)
HemeE16_regress$Group_Anno2[grep(pattern = "TCell", HemeE16_regress$Group_Anno2, fixed = TRUE)] <- "TCells"

### set the ident of the object with the group info
HemeE16_regress <- SetIdent(object = HemeE16_regress,
                            cells = rownames(HemeE16_regress@meta.data),
                            value = HemeE16_regress$Group_Anno2)

### get DE genes between two groups
de_result2 <- FindMarkers(object = HemeE16_regress,
                          test.use = "wilcox",
                          ident.1 = "TCells",
                          ident.2 = "Others",
                          logfc.threshold = 0.2,
                          min.pct = 0.2)

### order the DE reuslt
de_result2 <- de_result2[order(de_result2$p_val_adj),]

### data frame
de_result2 <- data.frame(Gene_Symbol=rownames(de_result2),
                         de_result2,
                          stringsAsFactors = FALSE, check.names = FALSE)

### now get the DE genes that are up-regulated in small T cells and up-regulated in the t cells (vs the other cells)
small_t_degs <- de_result$Gene_Symbol[intersect(which(de_result$p_val_adj < 0.0001),
                                                which(de_result$avg_log2FC < 0))]
t_degs <- de_result2$Gene_Symbol[intersect(which(de_result2$p_val_adj < 0.0001),
                                           which(de_result2$avg_log2FC > 0))]

### print some
writeLines(paste(intersect(small_t_degs, t_degs)))


#
### DE & pathway analysis between E16.5 vs E18.5 based on endothelial clusters
#

### updated stroma RDS file
Robj2_path="./data/Combined_Seurat_Obj.RDS"
Updated_Seurat_Obj <- readRDS(file = Robj2_path)

### create new column for the analysis
Updated_Seurat_Obj@meta.data$Dev_Anno <- paste0(Updated_Seurat_Obj@meta.data$Development, "_",
                                                Updated_Seurat_Obj@meta.data$Annotation)

### determine necessary variables
time_points <- c("E16", "E18", "P0", "ADULT")
HSPC_populations <- c("LTHSC", "STHSC", "MPP2", "MPP3", "MPP4")

### set the ident of the object with Tissue type
Updated_Seurat_Obj <- SetIdent(object = Updated_Seurat_Obj,
                               cells = rownames(Updated_Seurat_Obj@meta.data),
                               value = Updated_Seurat_Obj@meta.data$Tissue)

### check whether the orders are the same
print(identical(names(Idents(object = Updated_Seurat_Obj)), rownames(Updated_Seurat_Obj@meta.data)))

### run UMAP
Updated_Seurat_Obj <- FindVariableFeatures(Updated_Seurat_Obj)
Updated_Seurat_Obj <- ScaleData(Updated_Seurat_Obj)
Updated_Seurat_Obj <- RunUMAP(Updated_Seurat_Obj, dims = 1:5)

### since we are only interested in ECs, just take ECs
Updated_Seurat_Obj2 <- subset(Updated_Seurat_Obj,
                              cells = rownames(Updated_Seurat_Obj@meta.data)[which(Updated_Seurat_Obj$Annotation == "ECs")])

### annotate E16 & E18 ECs
Updated_Seurat_Obj2$Annotation2 <- paste0(Updated_Seurat_Obj2$Annotation, "_", Updated_Seurat_Obj2$Development)
Updated_Seurat_Obj2 <- SetIdent(object = Updated_Seurat_Obj2,
                                cells = rownames(Updated_Seurat_Obj2@meta.data),
                                value = Updated_Seurat_Obj2$Annotation2)

### perform DE analysis between E16 ECs vs E18 ECs
de_result <- FindMarkers(object = Updated_Seurat_Obj2,
                         test.use = "wilcox",
                         ident.1 = "ECs_E16",
                         ident.2 = "ECs_E18",
                         logfc.threshold = 0.2,
                         min.pct = 0.2)

### order the DE reuslt
de_result <- de_result[order(de_result$p_val_adj),]

### data frame
de_result <- data.frame(Gene_Symbol=rownames(de_result),
                        de_result,
                        stringsAsFactors = FALSE, check.names = FALSE)

### write out the DE result
write.xlsx2(de_result,
            file = paste0("./results/",
                          "ECs_E16_vs_E18.xlsx"),
            sheetName = "DE_Result",
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

### pathway analysis with up-regulated genes & FDR < 0.05
de_entrez_ids <- mapIds(org.Mm.eg.db,
                        rownames(de_result)[intersect(which(de_result$p_val_adj < 0.05),
                                                      which(de_result$avg_log2FC > 0))],
                        "ENTREZID", "SYMBOL")
de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]

### GO & KEGG
pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                        org = "mouse", database = "GO",
                                        title = paste0("ECs_E16_vs_E18_Upregulated_Pathways"),
                                        displayNum = 10, imgPrint = TRUE,
                                        dir = "./results/")
pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                          org = "mouse", database = "KEGG",
                                          title = paste0("ECs_E16_vs_E18_Upregulated_Pathways"),
                                          displayNum = 10, imgPrint = TRUE,
                                          dir = "./results/")
write.xlsx2(pathway_result_GO, file = paste0("./results/GO_ECs_E16_vs_E18_Upregulated_Pathways.xlsx"),
            row.names = FALSE, sheetName = paste0("GO_Results"))
write.xlsx2(pathway_result_KEGG, file = paste0("./results/KEGG_ECs_E16_vs_E18_Upregulated_Pathways.xlsx"),
            row.names = FALSE, sheetName = paste0("KEGG_Results"))

### pathway analysis with down-regulated genes & FDR < 0.05
de_entrez_ids <- mapIds(org.Mm.eg.db,
                        rownames(de_result)[intersect(which(de_result$p_val_adj < 0.05),
                                                      which(de_result$avg_log2FC < 0))],
                        "ENTREZID", "SYMBOL")
de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]

### GO & KEGG
pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                        org = "mouse", database = "GO",
                                        title = paste0("ECs_E16_vs_E18_Downregulated_Pathways"),
                                        displayNum = 10, imgPrint = TRUE,
                                        dir = "./results/")
pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                          org = "mouse", database = "KEGG",
                                          title = paste0("ECs_E16_vs_E18_Downregulated_Pathways"),
                                          displayNum = 10, imgPrint = TRUE,
                                          dir = "./results/")
write.xlsx2(pathway_result_GO, file = paste0("./results/GO_ECs_E16_vs_E18_Downregulated_Pathways.xlsx"),
            row.names = FALSE, sheetName = paste0("GO_Results"))
write.xlsx2(pathway_result_KEGG, file = paste0("./results/KEGG_ECs_E16_vs_E18_Downregulated_Pathways.xlsx"),
            row.names = FALSE, sheetName = paste0("KEGG_Results"))







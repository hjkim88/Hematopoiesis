###
#   File name : GO_Enrichment_Analysis.R
#   Author    : Hyunjin Kim
#   Date      : Aug 5, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : GO Analyses on the clusters within each R object.
#               I’d like to know how the transcriptional programs of these immunophenotypic populations
#               change across development.
#   
#   Instruction
#               1. Source("GO_Enrichment_Analysis.R")
#               2. Run the function "go_analysis_trent" - specify the input directory and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_GO_Enrichment_Analysis.R/GO_Enrichment_Analysis.R")
#               > go_analysis_trent(Robj_path="C:/Users/hkim8/SJ/HSPC Analyses/HSPC Subsets/",
#                                   outputDir="./results/")
###

go_analysis_trent <- function(Robj_path="C:/Users/hkim8/SJ/HSPC Analyses/HSPC Subsets/",
                              outputDir="./results/GO_Enrichment/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  options(java.parameters = "-Xmx10240m")
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
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
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
                theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("KEGG ", title))
              
              png(paste0(dir, "kegg_", title, "_CB.png"), width = 2000, height = 1000)
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
                theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("GO ", title))
              
              png(paste0(dir, "go_", title, "_CB.png"), width = 2000, height = 1000)
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
  
  ### get Robj file list
  f_list <- list.files(Robj_path, pattern = ".Robj$")
  
  ### for each of the Robj file, run the GO analysis
  for(f in f_list) {
    
    ### load the Seurat object and save the object name
    tmp_env <- new.env()
    load(paste0(Robj_path, f), tmp_env)
    obj_name <- ls(tmp_env)
    assign("Seurat_Obj", get(obj_name, envir = tmp_env))
    rm(tmp_env)
    gc()
    
    ### rownames in the meta.data should be in the same order as colnames in the counts
    Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
    
    ### Find markers for each cluster?
    ### Or FindAllMarkers()?
    
    ### print out the UMAP of the R object
    outputDir2 <- paste0(outputDir, f, "/")
    dir.create(outputDir2, recursive = TRUE, showWarnings = FALSE)
    if(nrow(Seurat_Obj@meta.data) > 1) {
      DimPlot(Seurat_Obj, reduction = "umap", group.by = "seurat_clusters", pt.size = 2) +
        labs(title = paste0("UMAP_", f, "_Clusters"))
      ggsave(file = paste0(outputDir2, "UMAP_", f, "_Clusters.png"), width = 12, height = 8, dpi = 300)
    }
    
    ### existing clusters
    cluster_list <- intersect(levels(Seurat_Obj@meta.data$seurat_clusters),
                              as.character(unique(Seurat_Obj@meta.data$seurat_clusters)))
    
    ### print the number of clusters in the R object
    writeLines(paste(f, "-Unique Cluster #-", length(cluster_list)))
    
    ### run only if there are more than one clusters
    if(length(cluster_list) > 1) {
      ### Ident configure
      Idents(object = Seurat_Obj) <- Seurat_Obj@meta.data$seurat_clusters
      
      ### DE analysis with DESeq2
      all_de_result <- FindAllMarkers(object = Seurat_Obj, test.use = "DESeq2",
                                      logfc.threshold = 0, min.pct = 0.1)
      
      ### GO Enrichment for each cluster
      for(cluster in cluster_list) {
        
        ### new output dir
        outputDir2 <- paste0(outputDir2, cluster, "/")
        dir.create(outputDir2, recursive = TRUE, showWarnings = FALSE)
        
        ### extract certain DE result
        de_result <- all_de_result[which(all_de_result$cluster == cluster),]
        
        if(nrow(de_result) > 0) {
          ### change the column order
          de_result <- data.frame(gene=de_result$gene,
                                  cluster=de_result$cluster,
                                  de_result[,-c(6,7)],
                                  stringsAsFactors = FALSE, check.names = FALSE)
          
          ### write out the result
          write.xlsx2(de_result, file = paste0(outputDir2, f, "_GO_enrichment_result_", cluster, ".xlsx"),
                      row.names = FALSE, sheetName = "DE_Result", append = FALSE)
          
          ### GO Enrichment analysis with the DE genes
          go_result <- pathwayAnalysis_CP(geneList = mapIds(org.Mm.eg.db, de_result$gene, "ENTREZID", "SYMBOL"),
                                          org = "mouse", database = "GO",
                                          title = paste0(f, "_results_", cluster),
                                          displayNum = 50, imgPrint = TRUE,
                                          dir = paste0(outputDir2))
          
          ### write out the GO result
          write.xlsx2(go_result, file = paste0(outputDir2, f, "_GO_enrichment_result_", cluster, ".xlsx"),
                      row.names = FALSE, sheetName = "GO_Result", append = TRUE)
        }
        
      }
    }
  }
  
}

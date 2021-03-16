###
#   File name : RNAMagnet_Pv_Calculation.R
#   Author    : Hyunjin Kim
#   Date      : Mar 15, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Make a function similar to RNAMagnetAnchors but also computes permutation p-values
#               based on randomly generated specificity scores
#   
#   Instruction
#               1. Source("RNAMagnet_Pv_Calculation.R")
#               2. Run the function "rnamagnet_pv_cal" - specify the input path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_RNAMagnet_Pv_Calculation.R/RNAMagnet_Pv_Calculation.R")
#               > rnamagnet_pv_cal(Robj_path="./data/New_Objects/AdultStromavLTHSC_NA.Robj",
#                                  outputDir="./results/New_Objects/")
###

rnamagnet_pv_cal <- function(Robj_path="./data/New_Objects/AdultStromavLTHSC_NA.Robj",
                             outputDir="./results/New_Objects/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(RNAMagnet, quietly = TRUE)) {
    remotes::install_github("veltenlab/rnamagnet")
    require(RNAMagnet, quietly = TRUE)
  }
  if(!require(reticulate, quietly = TRUE)) {
    install.packages("reticulate")
    require(reticulate, quietly = TRUE)
  }
  if(!require(Rmagic, quietly = TRUE)) {
    install.packages("Rmagic")
    require(Rmagic, quietly = TRUE)
  }
  use_condaenv("r-reticulate")
  py_config()
  py_module_available("magic")
  
  ### loads an RData file, and returns it
  loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
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
  RNAMagnetBase <- function(seurat, anchors=NULL,neighborhood.distance=NULL, neighborhood.gradient =NULL, .k = 10, .x0 = 0.5, .minExpression,.minMolecules=1, .version = "1.0.0", .cellularCompartment, .manualAnnotation = "Correct", .symmetric = F) {
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
  RNAMagnetAnchors <- function(seurat, anchors, return = "summary", neighborhood.distance = 0.7, neighborhood.gradient = 3, .k = 10, .x0 = 0.5, .minExpression = 0, .minMolecules = 1, .version = "1.0.0", .cellularCompartment = c("Membrane","ECM","Both"), .manualAnnotation = "Correct" ) {
    
    myMagnet <- RNAMagnetBase(seurat, anchors, neighborhood.distance,neighborhood.gradient, .k, .x0, .minExpression, .minMolecules, .version, .cellularCompartment, .manualAnnotation,TRUE)
    if (return=="rnamagnet-class") myMagnet else data.frame(direction = as.factor(colnames(myMagnet@specificity)[apply(myMagnet@specificity,1,which.max)]), adhesiveness = myMagnet@adhesiveness, myMagnet@specificity[,anchors])
    
  }
  
  ### RNAMagnetSignaling
  RNAMagnetSignaling <- function(seurat, neighborhood.distance = NULL, neighborhood.gradient = NULL, .k = 10, .x0 = 0.5, .minExpression = 10, .minMolecules = 1, .version = "1.0.0", .cellularCompartment = c("Secreted","Both"), .manualAnnotation = "Correct" ) {
    
    RNAMagnetBase(seurat, anchors = NULL, neighborhood.distance,neighborhood.gradient, .k, .x0, .minExpression, .minMolecules, .version, .cellularCompartment, .manualAnnotation, FALSE)
    
  }
  
  ### load the R object
  Seurat_Obj <- loadRData(Robj_path)
  gc()
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
  print(identical(rownames(Seurat_Obj@meta.data), colnames(Seurat_Obj@assays$RNA@counts)))
  
  ### active assay = "RNA"
  Seurat_Obj@active.assay <- "RNA"
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data)))
  
  ### set anchor column
  target_col <- "HSPC"
  
  ### RNAMagnet
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data[,target_col])
  RNAMagnet_Anchors_Result <- RNAMagnetAnchors(seurat = Seurat_Obj,
                                               anchors = unique(Seurat_Obj@meta.data[,target_col]),
                                               return = "summary",
                                               neighborhood.distance = 0.7,
                                               neighborhood.gradient = 3,
                                               .k = 10,
                                               .x0 = 0.5,
                                               .minExpression = 0,
                                               .minMolecules = 1,
                                               .version = "1.0.0",
                                               .cellularCompartment = c("Membrane","ECM","Both"),
                                               .manualAnnotation = "Correct")
  
  ### calculate specificity scores
  RNAMagnet_Anchors_Result$specificity <- as.numeric(sapply(rownames(RNAMagnet_Anchors_Result), function(x) {
    return(RNAMagnet_Anchors_Result[x, as.character(RNAMagnet_Anchors_Result[x,"direction"])])
  }))
  
  
  ### permutation p-values
  iteration_num <- 10
  set.seed(2990)
  random_value_list <- NULL
  unique_anchors <- unique(Seurat_Obj@meta.data[,target_col])
  unique_anchor_nums <- sapply(unique_anchors, function(x) length(which(Seurat_Obj@meta.data[,target_col] == x)))
  for(i in 1:iteration_num) {
    
    ### set random anchors similar to the original
    random_anchors <- rep("", sum(unique_anchor_nums))
    for(anchor in unique_anchors) {
      random_anchors[sample(which(random_anchors == ""), unique_anchor_nums[anchor])] <- anchor
    }
    
    ### run RNAMagnet signaling with random anchors
    Seurat_Obj <- SetIdent(object = Seurat_Obj,
                           cells = rownames(Seurat_Obj@meta.data),
                           value = random_anchors)
    random_result <- RNAMagnetAnchors(seurat = Seurat_Obj,
                                      anchors = unique(Seurat_Obj@meta.data[,target_col]),
                                      return = "summary",
                                      neighborhood.distance = 0.7,
                                      neighborhood.gradient = 3,
                                      .k = 10,
                                      .x0 = 0.5,
                                      .minExpression = 0,
                                      .minMolecules = 1,
                                      .version = "1.0.0",
                                      .cellularCompartment = c("Membrane","ECM","Both"),
                                      .manualAnnotation = "Correct")
    
    ### calculate specificity scores
    random_result$specificity <- as.numeric(sapply(rownames(random_result), function(x) {
      return(random_result[x, as.character(random_result[x,"direction"])])
    }))
    
    ### combine results from all the iterations
    if(is.null(random_value_list)) {
      random_value_list <- random_result$specificity
    } else {
      random_value_list <- c(random_value_list, random_result$specificity)
    }
    
  }
  random_value_list <- random_value_list[order(-random_value_list)]
  
  ### give permutation p-value
  RNAMagnet_Anchors_Result$Specificity_P_Val <- NA
  for(i in 1:nrow(RNAMagnet_Anchors_Result)) {
    RNAMagnet_Anchors_Result$Specificity_P_Val[i] <- round(length(which(random_value_list >= RNAMagnet_Anchors_Result$specificity[i]))/length(random_value_list),
                                                           digits = 8)
  }
  
  ### write out the result
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  write.xlsx2(data.frame(Cell=rownames(RNAMagnet_Anchors_Result), RNAMagnet_Anchors_Result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, paste0("RNAMagnet_Result_",
                                              substring(basename(Robj_path), 1, nchar(basename(Robj_path))-5),
                                              ".xlsx")),
              sheetName = paste("RNAMagnet_Result_",
                                paste(unique_anchors, collapse = "_")),
              row.names = FALSE)
  gc()
  
}

### LTHSC vs Adult Stroma
rnamagnet_pv_cal(Robj_path="./data/New_Objects/AdultStromavLTHSC_NA.Robj",
                 outputDir="./results/New_Objects/")

### MPP2 vs Adult Stroma
rnamagnet_pv_cal(Robj_path="./data/New_Objects/AdultStromavMPP2_NA.Robj",
                 outputDir="./results/New_Objects/")

### LTHSC vs E16 Stroma
rnamagnet_pv_cal(Robj_path="./data/New_Objects/E16StromavLTHSC_NA.Robj",
                 outputDir="./results/New_Objects/")

### MPP2 vs E16 Stroma
rnamagnet_pv_cal(Robj_path="./data/New_Objects/E16StromavMPP2_NA.Robj",
                 outputDir="./results/New_Objects/")

### LTHSC vs E18 Stroma
rnamagnet_pv_cal(Robj_path="./data/New_Objects/E18StromavLTHSC_NA.Robj",
                 outputDir="./results/New_Objects/")

### MPP2 vs E18 Stroma
rnamagnet_pv_cal(Robj_path="./data/New_Objects/E18StromavMPP2_NA.Robj",
                 outputDir="./results/New_Objects/")

### LTHSC vs P0 Stroma
rnamagnet_pv_cal(Robj_path="./data/New_Objects/P0StromavLTHSC_NA.Robj",
                 outputDir="./results/New_Objects/")

### MPP2 vs P0 Stroma
rnamagnet_pv_cal(Robj_path="./data/New_Objects/P0StromavMPP2_NA.Robj",
                 outputDir="./results/New_Objects/")

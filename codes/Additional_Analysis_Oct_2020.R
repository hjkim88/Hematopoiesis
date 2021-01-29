###
#   File name : Additional_Analysis_Oct_2020.R
#   Author    : Hyunjin Kim
#   Date      : Oct 28, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Perform additional analyses for Trent's project based on the
#               "20201012 Analyses for Hyunjin Overview.xlsx".
#   
#   Instruction
#               1. Source("Additional_Analysis_Oct_2020.R")
#               2. Run the function "additional_analysis" - specify the input paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Additional_Analysis_Oct_2020.R/Additional_Analysis_Oct_2020.R")
#               > additional_analysis(Robj1_path="./data/Combined_Seurat_Obj.RDATA",
#                                     Robj2_path="./data/Combined_Seurat_Obj.RDS",
#                                     outputDir="./results/Additional_Oct2020/")
###

additional_analysis <- function(Robj1_path="./data/Combined_Seurat_Obj.RDATA",
                                Robj2_path="./data/Combined_Seurat_Obj.RDS",
                                outputDir="./results/Additional_Oct2020/") {
  
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
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  if(!require(ggbeeswarm, quietly = TRUE)) {
    install.packages("ggbeeswarm")
    require(ggbeeswarm, quietly = TRUE)
  }
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
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
               cex = cex, pch = 16, col = col[1])
        if(show.constraints && !is.null(constraints.col)){
          for(const in names(constraints.col)) {
            points(centers[clusters %in% const, dims,
                           drop=FALSE], cex = cex / 2,
                   col = constraints.col[const], pch = 16)
            text(x = centers[clusters %in% const, dims[1]]+0,
                 y = centers[clusters %in% const, dims[2]]+8,
                 labels = const,
                 cex = cex / 2,
                 col = "black")
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
  
  #' @title Pairs plot of Slingshot output
  #' @name pairs-SlingshotDataSet
  #'
  #' @description A tool for quickly visualizing lineages inferred by
  #'   \code{slingshot}.
  #'
  #' @param x a \code{SlingshotDataSet} with results to be plotted.
  #' @param type character, the type of output to be plotted, can be one of
  #'   \code{"lineages"}, \code{curves}, or \code{both} (by partial matching), see
  #'   Details for more.
  #' @param show.constraints logical, whether or not the user-specified initial
  #'   and terminal clusters should be specially denoted by green and red dots,
  #'   respectively.
  #' @param col character, color vector for points.
  #' @param pch integer or character specifying the plotting symbol, see
  #'   \code{\link{par}}.
  #' @param cex numeric, amount by which points should be magnified, see
  #'   \code{\link{par}}.
  #' @param lwd numeric, the line width, see \code{\link{par}}.
  #' @param ... additional parameters for \code{plot} or \code{axis}, see
  #'   \code{\link{pairs}}.
  #' @param labels character, the names of the variables, see \code{\link{pairs}}.
  #' @param horInd see \code{\link{pairs}}.
  #' @param verInd see \code{\link{pairs}}.
  #' @param lower.panel see \code{\link{pairs}}.
  #' @param upper.panel see \code{\link{pairs}}.
  #' @param diag.panel see \code{\link{pairs}}.
  #' @param text.panel see \code{\link{pairs}}.
  #' @param label.pos see \code{\link{pairs}}.
  #' @param line.main see \code{\link{pairs}}.
  #' @param cex.labels see \code{\link{pairs}}.
  #' @param font.labels see \code{\link{pairs}}.
  #' @param row1attop see \code{\link{pairs}}.
  #' @param gap see \code{\link{pairs}}.
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
  #' pairs(sds, type = 'curves')
  #'
  #' @export
  pairs.SlingshotDataSet <-
    function (x, type = NULL, show.constraints = FALSE, col = NULL,
              constraints.col = NULL,
              pch = 16, cex=1, lwd=2, ...,
              labels, horInd = seq_len(nc), verInd = seq_len(nc),
              lower.panel = FALSE, upper.panel = TRUE,
              diag.panel = NULL, text.panel = textPanel,
              label.pos = 0.5 + has.diag/3, line.main = 3,
              cex.labels = NULL, font.labels = 1,
              row1attop = TRUE, gap = 1) {
      #####
      lp.sling <- lower.panel
      up.sling <- upper.panel
      panel <- points
      if(!up.sling){
        upper.panel <- NULL
      }else{
        upper.panel <- panel
      }
      if(!lower.panel){
        lower.panel <- NULL
      }else{
        lower.panel <- panel
      }
      log = ""
      sds <- x
      x <- reducedDim(sds)
      curves <- FALSE
      lineages <- FALSE
      if(is.null(type)){
        if(length(slingCurves(sds)) > 0){
          type <- 'curves'
        }else if(length(slingLineages(sds)) > 0){
          type <- 'lineages'
        }else{
          stop('No lineages or curves detected.')
        }
      }else{
        type <- c('curves','lineages','both')[pmatch(type,
                                                     c('curves','lineages',
                                                       'both'))]
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
      if(lineages & (length(slingLineages(sds))==0)){
        stop('No lineages detected.')
      }
      if(curves & (length(slingCurves(sds))==0)){
        stop('No curves detected.')
      }
      if(lineages){
        forest <- slingAdjacency(sds)
        clusters <- rownames(forest)
        nclus <- nrow(forest)
        centers <- t(vapply(clusters,function(clID){
          w <- slingClusterLabels(sds)[,clID]
          return(apply(x, 2, weighted.mean, w = w))
        }, rep(0,ncol(reducedDim(sds)))))
        rownames(centers) <- clusters
        linC <- slingParams(sds)
      }
      range.max <- max(apply(x,2,function(xi){
        r <- range(xi, na.rm = TRUE)
        return(abs(r[2] - r[1]))
      }))
      plot.ranges <- apply(x,2,function(xi){
        mid <- (max(xi,na.rm = TRUE) + min(xi,na.rm = TRUE))/2
        return(c(mid - range.max/2, mid + range.max/2))
      })
      if(is.null(col)){
        if(requireNamespace("RColorBrewer", quietly = TRUE)) {
          cc <- c(RColorBrewer::brewer.pal(9, "Set1")[-c(1,3,6)],
                  RColorBrewer::brewer.pal(7, "Set2")[-2],
                  RColorBrewer::brewer.pal(6, "Dark2")[-5],
                  RColorBrewer::brewer.pal(8, "Set3")[-c(1,2)])
        } else {
          cc <- seq_len(100)
        }
        col <- cc[apply(slingClusterLabels(sds),1,which.max)]
      }
      #####
      if(doText <- missing(text.panel) || is.function(text.panel))
        textPanel <-
        function(x = 0.5, y = 0.5, txt, cex, font)
          text(x, y, txt, cex = cex, font = font)
      
      localAxis <- function(side, x, y, xpd, bg, col=NULL, lwd=NULL, main,
                            oma, ...) {
        ## Explicitly ignore any color argument passed in as
        ## it was most likely meant for the data points and
        ## not for the axis.
        xpd <- NA
        if(side %% 2L == 1L && xl[j]) xpd <- FALSE
        if(side %% 2L == 0L && yl[i]) xpd <- FALSE
        if(side %% 2L == 1L) Axis(x, side = side, xpd = xpd, ...)
        else Axis(y, side = side, xpd = xpd, ...)
      }
      
      localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
      localLowerPanel <- function(..., main, oma, font.main, cex.main)
        lower.panel(...)
      localUpperPanel <- function(..., main, oma, font.main, cex.main)
        upper.panel(...)
      localDiagPanel <- function(..., main, oma, font.main, cex.main)
        diag.panel(...)
      
      dots <- list(...); nmdots <- names(dots)
      if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for(i in seq_along(names(x))) {
          if(is.factor(x[[i]]) || is.logical(x[[i]]))
            x[[i]] <- as.numeric(x[[i]])
          if(!is.numeric(unclass(x[[i]])))
            stop("non-numeric argument to 'pairs'")
        }
      } else if (!is.numeric(x)) stop("non-numeric argument to 'pairs'")
      panel <- match.fun(panel)
      if((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
        lower.panel <- match.fun(lower.panel)
      if((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
        upper.panel <- match.fun(upper.panel)
      if((has.diag  <- !is.null( diag.panel)) && !missing( diag.panel))
        diag.panel <- match.fun( diag.panel)
      
      if(row1attop) {
        tmp <- lower.panel; lower.panel <- upper.panel; upper.panel <- tmp
        tmp <- has.lower; has.lower <- has.upper; has.upper <- tmp
      }
      
      nc <- ncol(x)
      if (nc < 2L) stop("only one column in the argument to 'pairs'")
      if(!all(horInd >= 1L & horInd <= nc))
        stop("invalid argument 'horInd'")
      if(!all(verInd >= 1L & verInd <= nc))
        stop("invalid argument 'verInd'")
      if(doText) {
        if (missing(labels)) {
          labels <- colnames(x)
          if (is.null(labels)) labels <- paste("var", 1L:nc)
        }
        else if(is.null(labels)) doText <- FALSE
      }
      oma <- if("oma" %in% nmdots) dots$oma
      main <- if("main" %in% nmdots) dots$main
      if (is.null(oma))
        oma <- c(4, 4, if(!is.null(main)) 6 else 4, 4)
      opar <- par(mfrow = c(length(horInd), length(verInd)),
                  mar = rep.int(gap/2, 4), oma = oma)
      on.exit(par(opar))
      dev.hold(); on.exit(dev.flush(), add = TRUE)
      
      xl <- yl <- logical(nc)
      if (is.numeric(log)) xl[log] <- yl[log] <- TRUE
      else {xl[] <- grepl("x", log); yl[] <- grepl("y", log)}
      for (i in if(row1attop) verInd else rev(verInd))
        for (j in horInd) {
          l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", ""))
          localPlot(x[, j], x[, i], xlab = "", ylab = "",
                    axes = FALSE, type = "n", ..., log = l,
                    xlim = plot.ranges[,j], ylim = plot.ranges[,i])
          if(i == j || (i < j && has.lower) || (i > j && has.upper) ) {
            box()
            if(i == 1  && (!(j %% 2L) || !has.upper || !has.lower ))
              localAxis(1L + 2L*row1attop, x[, j], x[, i], ...)
            if(i == nc && (  j %% 2L  || !has.upper || !has.lower ))
              localAxis(3L - 2L*row1attop, x[, j], x[, i], ...)
            if(j == 1  && (!(i %% 2L) || !has.upper || !has.lower ))
              localAxis(2L, x[, j], x[, i], ...)
            if(j == nc && (  i %% 2L  || !has.upper || !has.lower ))
              localAxis(4L, x[, j], x[, i], ...)
            mfg <- par("mfg")
            if(i == j) {
              if (has.diag) localDiagPanel(as.vector(x[, i]), ...)
              if (doText) {
                par(usr = c(0, 1, 0, 1))
                if(is.null(cex.labels)) {
                  l.wid <- strwidth(labels, "user")
                  cex.labels <- max(0.8, min(2, .9 / max(l.wid)))
                }
                xlp <- if(xl[i]) 10^0.5 else 0.5
                ylp <- if(yl[j]) 10^label.pos else label.pos
                text.panel(xlp, ylp, labels[i],
                           cex = cex.labels, font = font.labels)
              }
            } else if(i < j){
              if(up.sling){
                points(as.vector(x[, j]), as.vector(x[, i]),
                       col = col, cex = cex, pch=pch, ...)
                if(lineages){
                  for(ii in seq_len(nclus-1)){
                    for(jj in seq(ii+1,nclus)){
                      if(forest[ii,jj]==1){
                        seg.col <- 1
                        lines(centers[c(ii,jj),j],
                              centers[c(ii,jj),i],
                              lwd = lwd, col = seg.col, ...)
                      }
                    }
                  }
                  points(centers[,j],centers[,i], pch = pch,
                         cex=2*cex)
                  if(show.constraints && is.null(constraints.col)){
                    if(any(linC$start.given)){
                      st.ind <- clusters %in%
                        linC$start.clus[linC$start.given]
                      points(centers[st.ind,j],
                             centers[st.ind,i], cex = cex,
                             col = 'green3',
                             pch = pch)
                    }
                    if(any(linC$end.given)){
                      en.ind <- clusters %in%
                        linC$end.clus[linC$end.given]
                      points(centers[en.ind,j],
                             centers[en.ind,i], cex = cex,
                             col = 'red2', pch = pch)
                    }
                  } else if(show.constraints && !is.null(constraints.col)){
                    for(const in names(constraints.col)) {
                      points(centers[clusters %in% const, j, drop=FALSE],
                             centers[clusters %in% const, i, drop=FALSE],
                             cex = cex, pch = 16,
                             col = constraints.col[const])
                    }
                  }
                }
                if(curves){
                  for(c in slingCurves(sds)){
                    lines(c$s[c$ord,c(j,i)], lwd = lwd,
                          col=1, ...)
                  }
                }
              }
            }
            else{
              if(lp.sling){
                points(as.vector(x[, j]), as.vector(x[, i]),
                       col = col, cex = cex, pch=pch, ...)
                if(lineages){
                  for(ii in seq_len(nclus-1)){
                    for(jj in seq(ii+1,nclus)){
                      if(forest[ii,jj]==1){
                        if(clusters[ii] %in%
                           linC$start.clus |
                           clusters[jj] %in%
                           linC$start.clus){
                          seg.col <- 'green3'
                        }else if(clusters[ii] %in%
                                 linC$end.clus[
                                   linC$end.given] |
                                 clusters[jj] %in%
                                 linC$end.clus[
                                   linC$end.given]){
                          seg.col <- 'red2'
                        }else{
                          seg.col <- 1
                        }
                        lines(centers[c(ii,jj),j],
                              centers[c(ii,jj),i],
                              lwd = lwd, col = seg.col,...)
                      }
                    }
                  }
                  points(centers[,j],centers[,i], pch = pch,
                         cex = 2*cex)
                }
                if(curves){
                  for(c in slingCurves(sds)){
                    lines(c$s[c$ord,c(j,i)],lwd = lwd,
                          col=1, ...)
                  }
                }
              }
            }
            if (any(par("mfg") != mfg))
              stop("the 'panel' function made a new plot")
          } else par(new = FALSE)
          
        }
      if (!is.null(main)) {
        font.main <- if("font.main" %in% nmdots){
          dots$font.main
        }else par("font.main")
        cex.main <- if("cex.main" %in%
                       nmdots) dots$cex.main else par("cex.main")
        mtext(main, 3, line.main, outer=TRUE, at = 0.5, cex = cex.main,
              font = font.main)
      }
      invisible(NULL)
  }
  
  
  ###
  #   This function was downloaded from: https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R
  ###
  heatmap.3 <- function(x,
                        Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                        distfun = dist,
                        hclustfun = hclust,
                        dendrogram = c("both","row", "column", "none"),
                        symm = FALSE,
                        scale = c("none","row", "column"),
                        na.rm = TRUE,
                        revC = identical(Colv,"Rowv"),
                        add.expr,
                        breaks,
                        symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                        col = "heat.colors",
                        colsep,
                        rowsep,
                        sepcolor = "white",
                        sepwidth = c(0.05, 0.05),
                        cellnote,
                        notecex = 1,
                        notecol = "cyan",
                        na.color = par("bg"),
                        trace = c("none", "column","row", "both"),
                        tracecol = "cyan",
                        hline = median(breaks),
                        vline = median(breaks),
                        linecol = tracecol,
                        margins = c(5,5),
                        ColSideColors,
                        RowSideColors,
                        side.height.fraction=0.3,
                        cexRow = 0.2 + 1/log10(nr),
                        cexCol = 0.2 + 1/log10(nc),
                        labRow = NULL,
                        labCol = NULL,
                        key = TRUE,
                        keysize = 1.5,
                        density.info = c("none", "histogram", "density"),
                        denscol = tracecol,
                        symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                        densadj = 0.25,
                        main = NULL,
                        xlab = NULL,
                        ylab = NULL,
                        lmat = NULL,
                        lhei = NULL,
                        lwid = NULL,
                        ColSideColorsSize = 1,
                        RowSideColorsSize = 1,
                        KeyValueName="Value",...){
    
    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
        return(TRUE)
      if (is.list(x))
        return(all(sapply(x, invalid)))
      else if (is.vector(x))
        return(all(is.na(x)))
      else return(FALSE)
    }
    
    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
      x <- (x - low)/(high - low)
      x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
      "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
      col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
      warning("Using scale=\"row\" or scale=\"column\" when breaks are",
              "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
      Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
      Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
      Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
      stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
      stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
      stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
      cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
      if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                   c("both", "row"))) {
        if (is.logical(Colv) && (Colv))
          dendrogram <- "column"
        else dedrogram <- "none"
        warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
      }
    }
    if (!inherits(Colv, "dendrogram")) {
      if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                   c("both", "column"))) {
        if (is.logical(Rowv) && (Rowv))
          dendrogram <- "row"
        else dendrogram <- "none"
        warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
      }
    }
    if (inherits(Rowv, "dendrogram")) {
      ddr <- Rowv
      rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
      Rowv <- rowMeans(x, na.rm = na.rm)
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else {
      rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
      ddc <- Colv
      colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      if (exists("ddr")) {
        ddc <- ddr
        colInd <- order.dendrogram(ddc)
      }
      else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
      Colv <- colMeans(x, na.rm = na.rm)
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else {
      colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
      labRow <- if (is.null(rownames(x)))
        (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
      labCol <- if (is.null(colnames(x)))
        (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
      retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
      x <- sweep(x, 1, rm)
      retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
      x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
      retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
      x <- sweep(x, 2, rm)
      retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
      x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
      if (missing(col) || is.function(col))
        breaks <- 16
      else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
      if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                      length = breaks)
      else {
        extreme <- max(abs(x), na.rm = TRUE)
        breaks <- seq(-extreme, extreme, length = breaks)
      }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
      col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
      lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
      lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
      lmat <- rbind(4:3, 2:1)
      
      if (!missing(ColSideColors)) {
        #if (!is.matrix(ColSideColors))
        #stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
          stop("'ColSideColors' must be a matrix of nrow(x) rows")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        #lhei <- c(lhei[1], 0.2, lhei[2])
        lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
      }
      
      if (!missing(RowSideColors)) {
        #if (!is.matrix(RowSideColors))
        #stop("'RowSideColors' must be a matrix")
        if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
          stop("'RowSideColors' must be a matrix of ncol(x) columns")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
        #lwid <- c(lwid[1], 0.2, lwid[2])
        lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
      }
      lmat[is.na(lmat)] <- 0
    }
    
    if (length(lhei) != nrow(lmat))
      stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
      stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    if (!missing(RowSideColors)) {
      if (!is.matrix(RowSideColors)){
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
      } else {
        par(mar = c(margins[1], 0, 0, 0.5))
        rsc = t(RowSideColors[,rowInd, drop=F])
        rsc.colors = matrix()
        rsc.names = names(table(rsc))
        rsc.i = 1
        for (rsc.name in rsc.names) {
          rsc.colors[rsc.i] = rsc.name
          rsc[rsc == rsc.name] = rsc.i
          rsc.i = rsc.i + 1
        }
        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (length(rownames(RowSideColors)) > 0) {
          axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE, cex.axis = cexRow/2)
        }
      }
    }
    
    if (!missing(ColSideColors)) {
      
      if (!is.matrix(ColSideColors)){
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
      } else {
        par(mar = c(0.5, 0, 0, margins[2]))
        csc = ColSideColors[colInd, , drop=F]
        csc.colors = matrix()
        csc.names = names(table(csc))
        csc.i = 1
        for (csc.name in csc.names) {
          csc.colors[csc.i] = csc.name
          csc[csc == csc.name] = csc.i
          csc.i = csc.i + 1
        }
        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
        image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (length(colnames(ColSideColors)) > 0) {
          axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE, cex.axis = cexCol/2)
        }
      }
    }
    
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
      iy <- nr:1
      if (exists("ddr"))
        ddr <- rev(ddr)
      x <- x[, iy]
      cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
      retval$rowDendrogram <- ddr
    if (exists("ddc"))
      retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
      mmat <- ifelse(is.na(x), 1, NA)
      image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
         cex.axis = cexCol)
    if (!is.null(xlab))
      mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
         cex.axis = cexRow)
    if (!is.null(ylab))
      mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
      eval(substitute(add.expr))
    if (!missing(colsep))
      for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
      for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
      retval$vline <- vline
      vline.vals <- scale01(vline, min.scale, max.scale)
      for (i in colInd) {
        if (!is.null(vline)) {
          abline(v = i - 0.5 + vline.vals, col = linecol,
                 lty = 2)
        }
        xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
        xv <- c(xv[1], xv)
        yv <- 1:length(xv) - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (trace %in% c("both", "row")) {
      retval$hline <- hline
      hline.vals <- scale01(hline, min.scale, max.scale)
      for (i in rowInd) {
        if (!is.null(hline)) {
          abline(h = i + hline, col = linecol, lty = 2)
        }
        yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
        yv <- rev(c(yv[1], yv))
        xv <- length(yv):1 - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (!missing(cellnote))
      text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
           col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
      plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
      plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
      title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
      par(mar = c(5, 4, 2, 1), cex = 0.75)
      tmpbreaks <- breaks
      if (symkey) {
        max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
        min.raw <- -max.raw
        tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
        tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
      }
      else {
        min.raw <- min(x, na.rm = TRUE)
        max.raw <- max(x, na.rm = TRUE)
      }
      
      z <- seq(min.raw, max.raw, length = length(col))
      image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
      par(usr = c(0, 1, 0, 1))
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      axis(1, at = xv, labels = lv)
      if (scale == "row")
        mtext(side = 1, "Row Z-Score", line = 2)
      else if (scale == "column")
        mtext(side = 1, "Column Z-Score", line = 2)
      else mtext(side = 1, KeyValueName, line = 2)
      if (density.info == "density") {
        dens <- density(x, adjust = densadj, na.rm = TRUE)
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[-omit]
        dens$y <- dens$y[-omit]
        dens$x <- scale01(dens$x, min.raw, max.raw)
        lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
              lwd = 1)
        axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
        title("Color Key\nand Density Plot")
        par(cex = 0.5)
        mtext(side = 2, "Density", line = 2)
      }
      else if (density.info == "histogram") {
        h <- hist(x, plot = FALSE, breaks = breaks)
        hx <- scale01(breaks, min.raw, max.raw)
        hy <- c(h$counts, h$counts[length(h$counts)])
        lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
              col = denscol)
        axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
        title("Color Key\nand Histogram")
        par(cex = 0.5)
        mtext(side = 2, "Count", line = 2)
      }
      else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                    high = retval$breaks[-1], color = retval$col)
    invisible(retval)
  }
  
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
  
  
  ### hierarchical clustering functions
  dist.spear <- function(x) as.dist(1-cor(t(x), method = "spearman"))
  hclust.ave <- function(x) hclust(x, method="average")
  
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
                ggtitle(paste0("KEGG ", title)) +
                theme(axis.text = element_text(size = 25))
              
              png(paste0(dir, "kegg_", title, ".png"), width = 2000, height = 1000)
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
                ggtitle(paste0("GO ", title)) +
                theme(axis.text = element_text(size = 25))
              
              png(paste0(dir, "go_", title, ".png"), width = 2000, height = 1000)
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
  
  ### a function to get important & garbage genes, make heamtaps with those genes,
  ### and pathway analysis with those genes, GSEA with PC contribution,
  ### heatmaps with DE genes from comparison, pathway analysis with DE genes from
  ### comparison, GSEA with DE genes from comparison.
  ### Seurat_object: The seurat object to be analysed
  ### target ident: the ident that will be analyzed, if NULL, just use the given object without split
  ### target_col: a column name of the meta data of the seurat object that will be used in
  ###             the pseudotime analysis and in creating heatmaps. Usually a time-associated column
  ### target_col_factor_level: a factor level of the 'target_col'
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
    
    #
    ### heatmap
    #
    
    ### important genes
    
    ### set colside colors
    uniqueV <- unique(local_Seurat_Obj@meta.data[,target_col])
    colors <- colorRampPalette(brewer.pal(9,"Blues"))(length(uniqueV))
    names(colors) <- uniqueV
    
    ### get a matrix for the heatmap
    if(length(important_genes) > 300) {
      heatmap_mat <- data.frame(local_Seurat_Obj@assays$RNA@counts[important_genes[1:300],], check.names = FALSE)
    } else {
      heatmap_mat <- data.frame(local_Seurat_Obj@assays$RNA@counts[important_genes,], check.names = FALSE)
    }
    
    ### scale the data
    heatmap_mat_scaled <- scale_h(heatmap_mat, type = "row")
    
    ### because there are some outliers in positive values
    ### we set the maximum as abs(minimum)
    heatmap_mat_scaled[which(heatmap_mat_scaled > abs(min(heatmap_mat_scaled)))] <- abs(min(heatmap_mat_scaled))
    
    ### heatmap
    png(paste0(result_dir, PC, "_Heatmap_with_important_genes_", important_thresh, ".png"), width = 12000, height = 6000)
    par(oma=c(0,0,10,6))
    heatmap.3(as.matrix(heatmap_mat_scaled),
              xlab = "", ylab = "", col=greenred(300),
              scale="none", key=T, keysize=0.8, density.info="density",
              dendrogram = "none", trace = "none",
              labRow = rownames(heatmap_mat_scaled), labCol = FALSE,
              Rowv = TRUE, Colv = FALSE,
              distfun=dist.spear, hclustfun=hclust.ave,
              ColSideColors = cbind(colors[as.character(local_Seurat_Obj@meta.data[,target_col])]),
              cexRow = 1.9, cexCol = 1.9, na.rm = TRUE)
    title(main = paste0(PC, "_Genes_Heatmap_(",
                        nrow(heatmap_mat_scaled), " Genes x ",
                        ncol(heatmap_mat_scaled), " Cells)"),
          cex.main = 15, line = -10, outer = TRUE)
    legend("left", inset = 0, xpd = TRUE, title = "Time", legend = names(colors), fill = colors, cex = 15, box.lty = 0)
    dev.off()
    
    ### garbage genes
    
    ### get a matrix for the heatmap
    if(length(garbage_genes) > 300) {
      heatmap_mat <- data.frame(local_Seurat_Obj@assays$RNA@counts[garbage_genes[1:300],], check.names = FALSE)
    } else {
      heatmap_mat <- data.frame(local_Seurat_Obj@assays$RNA@counts[garbage_genes,], check.names = FALSE)
    }
    
    ### scale the data
    heatmap_mat_scaled <- scale_h(heatmap_mat, type = "row")
    
    ### because there are some outliers in positive values
    ### we set the maximum as abs(minimum)
    heatmap_mat_scaled[which(heatmap_mat_scaled > abs(min(heatmap_mat_scaled)))] <- abs(min(heatmap_mat_scaled))
    
    ### heatmap
    png(paste0(result_dir, PC, "_Heatmap_with_garbage_genes_",  garbage_thresh, ".png"), width = 12000, height = 6000)
    par(oma=c(0,0,10,6))
    heatmap.3(as.matrix(heatmap_mat_scaled),
              xlab = "", ylab = "", col=greenred(300),
              scale="none", key=T, keysize=0.8, density.info="density",
              dendrogram = "none", trace = "none",
              labRow = rownames(heatmap_mat_scaled), labCol = FALSE,
              Rowv = TRUE, Colv = FALSE,
              distfun=dist.spear, hclustfun=hclust.ave,
              ColSideColors = cbind(colors[as.character(local_Seurat_Obj@meta.data[,target_col])]),
              cexRow = 1.9, cexCol = 1.9, na.rm = TRUE)
    title(main = paste0(PC, "_Genes_Heatmap_(",
                        nrow(heatmap_mat_scaled), " Genes x ",
                        ncol(heatmap_mat_scaled), " Cells)"),
          cex.main = 15, line = -10, outer = TRUE)
    legend("left", inset = 0, xpd = TRUE, title = "Time", legend = names(colors), fill = colors, cex = 15, box.lty = 0)
    dev.off()
    
    
    ### pathway analysis with the important genes of the given PC
    if(species[1] == "human") {
      pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                                important_genes,
                                                                "ENTREZID", "SYMBOL"),
                                              org = species[1], database = "GO",
                                              title = paste0(PC, "_Pathway_Results_with_", length(important_genes), "_important_genes_", important_thresh),
                                              displayNum = 50, imgPrint = TRUE,
                                              dir = paste0(result_dir))
      pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                                  important_genes,
                                                                  "ENTREZID", "SYMBOL"),
                                                org = species[1], database = "KEGG",
                                                title = paste0(PC, "_Pathway_Results_with_", length(important_genes), "_important_genes_", important_thresh),
                                                displayNum = 50, imgPrint = TRUE,
                                                dir = paste0(result_dir))
    } else if(species[1] == "mouse") {
      pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Mm.eg.db,
                                                                important_genes,
                                                                "ENTREZID", "SYMBOL"),
                                              org = species[1], database = "GO",
                                              title = paste0(PC, "_Pathway_Results_with_", length(important_genes), "_important_genes_", important_thresh),
                                              displayNum = 50, imgPrint = TRUE,
                                              dir = paste0(result_dir))
      pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Mm.eg.db,
                                                                  important_genes,
                                                                  "ENTREZID", "SYMBOL"),
                                                org = species[1], database = "KEGG",
                                                title = paste0(PC, "_Pathway_Results_with_", length(important_genes), "_important_genes_", important_thresh),
                                                displayNum = 50, imgPrint = TRUE,
                                                dir = paste0(result_dir))
    }
    write.xlsx2(pathway_result_GO, file = paste0(result_dir, PC, "_GO_Pathway_Results_with_", length(important_genes), "_important_genes_", important_thresh, ".xlsx"),
                row.names = FALSE, sheetName = paste0("GO_Results"))
    write.xlsx2(pathway_result_KEGG, file = paste0(result_dir, PC, "_KEGG_Pathway_Results_with_", length(important_genes), "_important_genes_", important_thresh, ".xlsx"),
                row.names = FALSE, sheetName = paste0("KEGG_Results"))
    
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
                                              which(GSEA_result$NES > 2))]
    if(length(pathways) < 5) {
      pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-03),
                                                          which(GSEA_result$size > 30)),
                                                which(GSEA_result$NES > 1.8))]
    }
    if(length(pathways) < 5) {
      pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-03),
                                                          which(GSEA_result$size > 30)),
                                                which(GSEA_result$NES > 1.6))]
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
                                                       which(de_result$avg_logFC > 0))]
          down_de_genes <- rownames(de_result)[intersect(which(de_result$p_val_adj < 0.01),
                                                         which(de_result$avg_logFC < 0))]
          
          ### heatmap
          ### set colside colors
          uniqueV <- unique(local_Seurat_Obj@meta.data[,target_col])
          colors <- colorRampPalette(brewer.pal(9,"Blues"))(length(uniqueV)+1)[-1]
          names(colors) <- uniqueV
          
          ### get a matrix for the heatmap
          if(length(de_genes) > 300) {
            heatmap_mat <- data.frame(local_Seurat_Obj@assays$RNA@counts[de_genes[1:300],], check.names = FALSE)
          } else {
            heatmap_mat <- data.frame(local_Seurat_Obj@assays$RNA@counts[de_genes,], check.names = FALSE)
          }
          
          ### scale the data
          heatmap_mat_scaled <- scale_h(heatmap_mat, type = "row")
          
          ### because there are some outliers in positive values
          ### we set the maximum as abs(minimum)
          heatmap_mat_scaled[which(heatmap_mat_scaled > abs(min(heatmap_mat_scaled)))] <- abs(min(heatmap_mat_scaled))
          
          ### heatmap
          png(paste0(result_dir2, PC, "_", PC_Val, "_Comparison_Heatmap.png"), width = 12000, height = 6000)
          par(oma=c(0,0,10,6))
          heatmap.3(as.matrix(heatmap_mat_scaled),
                    xlab = "", ylab = "", col=greenred(300),
                    scale="none", key=T, keysize=0.8, density.info="density",
                    dendrogram = "none", trace = "none",
                    labRow = rownames(heatmap_mat_scaled), labCol = FALSE,
                    Rowv = TRUE, Colv = FALSE,
                    distfun=dist.spear, hclustfun=hclust.ave,
                    ColSideColors = cbind(colors[as.character(local_Seurat_Obj@meta.data[,target_col])]),
                    cexRow = 1.9, cexCol = 1.9, na.rm = TRUE)
          title(main = paste0(PC, "_Genes_Heatmap_(",
                              nrow(heatmap_mat_scaled), " Genes x ",
                              ncol(heatmap_mat_scaled), " Cells)"),
                cex.main = 15, line = -10, outer = TRUE)
          legend("left", inset = 0, xpd = TRUE, title = "Time", legend = names(colors), fill = colors, cex = 15, box.lty = 1, box.lwd = 5)
          dev.off()
          
          ### separate the heatmap as grp1 vs grp2
          heatmap_mat_scaled <- heatmap_mat_scaled[,c(grp1_idx, grp2_idx)]
          
          ### heatmap
          png(paste0(result_dir2, PC, "_", PC_Val, "_Comparison_Heatmap_Grouped.png"), width = 12000, height = 6000)
          par(oma=c(0,0,10,6))
          heatmap.3(as.matrix(heatmap_mat_scaled),
                    xlab = "", ylab = "", col=greenred(300),
                    scale="none", key=T, keysize=0.8, density.info="density",
                    dendrogram = "none", trace = "none",
                    labRow = rownames(heatmap_mat_scaled), labCol = FALSE,
                    Rowv = FALSE, Colv = FALSE,
                    distfun=dist.spear, hclustfun=hclust.ave,
                    ColSideColors = cbind(colors[as.character(local_Seurat_Obj@meta.data[c(grp1_idx, grp2_idx),target_col])],
                                          c(rep("cyan3", length(grp1_idx)), rep("maroon", length(grp2_idx)))),
                    cexRow = 1.9, cexCol = 1.9, na.rm = TRUE)
          title(main = paste0(PC, "_Genes_Heatmap_(",
                              nrow(heatmap_mat_scaled), " Genes x ",
                              ncol(heatmap_mat_scaled), " Cells)"),
                cex.main = 15, line = -10, outer = TRUE)
          legend("topright", inset = 0, xpd = TRUE, title = "Group", legend = c(paste0(PC, " >= ", PC_Val),
                                                                               paste0(PC, " < ", PC_Val)),
                 fill = c("cyan3", "maroon"), cex = 10, box.lty = 1, box.lwd = 5)
          legend("left", inset = 0, xpd = TRUE, title = "Time", legend = names(colors), fill = colors, cex = 15, box.lty = 1, box.lwd = 5)
          dev.off()
          
          ### with 50 +/- logFC genes
          grp1_num <- ifelse(length(up_de_genes) >= 50, 50, length(up_de_genes))
          grp2_num <- ifelse(length(down_de_genes) >= 50, 50, length(down_de_genes))
          heatmap_mat <- data.frame(local_Seurat_Obj@assays$RNA@counts[c(up_de_genes[1:grp1_num],
                                                                         down_de_genes[1:grp2_num]),],
                                    check.names = FALSE)
          heatmap_mat_scaled <- scale_h(heatmap_mat, type = "row")
          heatmap_mat_scaled[which(heatmap_mat_scaled > abs(min(heatmap_mat_scaled)))] <- abs(min(heatmap_mat_scaled))
          heatmap_mat_scaled <- heatmap_mat_scaled[,c(grp1_idx, grp2_idx)]
          
          ### heatmap
          png(paste0(result_dir2, PC, "_", PC_Val, "_Comparison_Heatmap_Grouped_50_logFCs.png"), width = 12000, height = 6000)
          par(oma=c(0,0,10,6))
          heatmap.3(as.matrix(heatmap_mat_scaled),
                    xlab = "", ylab = "", col=greenred(300),
                    scale="none", key=T, keysize=0.8, density.info="density",
                    dendrogram = "none", trace = "none",
                    labRow = rownames(heatmap_mat_scaled), labCol = FALSE,
                    Rowv = FALSE, Colv = FALSE,
                    distfun=dist.spear, hclustfun=hclust.ave,
                    ColSideColors = cbind(colors[as.character(local_Seurat_Obj@meta.data[c(grp1_idx, grp2_idx),target_col])],
                                          c(rep("cyan3", length(grp1_idx)), rep("maroon", length(grp2_idx)))),
                    cexRow = 1.9, cexCol = 1.9, na.rm = TRUE)
          title(main = paste0(PC, "_Genes_Heatmap_(",
                              nrow(heatmap_mat_scaled), " Genes x ",
                              ncol(heatmap_mat_scaled), " Cells)"),
                cex.main = 15, line = -10, outer = TRUE)
          legend("topright", inset = 0, xpd = TRUE, title = "Group", legend = c(paste0(PC, " >= ", PC_Val),
                                                                                paste0(PC, " < ", PC_Val)),
                 fill = c("cyan3", "maroon"), cex = 10, box.lty = 1, box.lwd = 5)
          legend("left", inset = 0, xpd = TRUE, title = "Time", legend = names(colors), fill = colors, cex = 15, box.lty = 1, box.lwd = 5)
          dev.off()
          
          
          ### pathway analysis
          if(species[1] == "human") {
            pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                                      de_genes,
                                                                      "ENTREZID", "SYMBOL"),
                                                    org = species[1], database = "GO",
                                                    title = paste0(PC, "_", PC_Val, "_Pathway_Results_with_", length(de_genes), "_DE_genes"),
                                                    displayNum = 50, imgPrint = TRUE,
                                                    dir = paste0(result_dir2))
            pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                                        de_genes,
                                                                        "ENTREZID", "SYMBOL"),
                                                      org = species[1], database = "KEGG",
                                                      title = paste0(PC, "_", PC_Val, "_Pathway_Results_with_", length(de_genes), "_DE_genes"),
                                                      displayNum = 50, imgPrint = TRUE,
                                                      dir = paste0(result_dir2))
          } else if(species[1] == "mouse") {
            pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Mm.eg.db,
                                                                      de_genes,
                                                                      "ENTREZID", "SYMBOL"),
                                                    org = species[1], database = "GO",
                                                    title = paste0(PC, "_", PC_Val, "_Pathway_Results_with_", length(de_genes), "_DE_genes"),
                                                    displayNum = 50, imgPrint = TRUE,
                                                    dir = paste0(result_dir2))
            pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Mm.eg.db,
                                                                        de_genes,
                                                                        "ENTREZID", "SYMBOL"),
                                                      org = species[1], database = "KEGG",
                                                      title = paste0(PC, "_", PC_Val, "_Pathway_Results_with_", length(de_genes), "_DE_genes"),
                                                      displayNum = 50, imgPrint = TRUE,
                                                      dir = paste0(result_dir2))
          }
          write.xlsx2(pathway_result_GO, file = paste0(result_dir2, PC, "_", PC_Val, "_GO_Pathway_Results_with_", length(de_genes), "_DE_genes.xlsx"),
                      row.names = FALSE, sheetName = paste0("GO_Results"))
          write.xlsx2(pathway_result_KEGG, file = paste0(result_dir2, PC, "_", PC_Val, "_KEGG_Pathway_Results_with_", length(de_genes), "_DE_genes.xlsx"),
                      row.names = FALSE, sheetName = paste0("KEGG_Results"))
          
          ### GSEA
          ### signature preparation
          signat <- de_result$avg_logFC
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
                                                      which(GSEA_result$NES > 1.8))]
          }
          if(length(pathways) < 5) {
            pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-03),
                                                                which(GSEA_result$size > 30)),
                                                      which(GSEA_result$NES > 1.6))]
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
  
  ### ALL
  
  ### set the ident of the object with the HSPC type
  Combined_Seurat_Obj <- SetIdent(object = Combined_Seurat_Obj,
                                  cells = rownames(Combined_Seurat_Obj@meta.data),
                                  value = Combined_Seurat_Obj@meta.data$HSPC)
  
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Combined_Seurat_Obj)), rownames(Combined_Seurat_Obj@meta.data)))
  
  ### new output directory
  type <- "MPP2"
  outputDir2 <- paste0(outputDir, "ALL/", type, "/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### split the Seurat obj based on HSPC info
  subset_Seurat_Obj <- subset(Combined_Seurat_Obj, idents=type)
  
  ### order the meta data by developmental time
  subset_Seurat_Obj@meta.data <- subset_Seurat_Obj@meta.data[order(factor(subset_Seurat_Obj@meta.data$Development,
                                                                          levels = time_points)),]
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  subset_Seurat_Obj@assays$RNA@counts <- subset_Seurat_Obj@assays$RNA@counts[,rownames(subset_Seurat_Obj@meta.data)]
  
  ### run PCA
  subset_Seurat_Obj <- RunPCA(subset_Seurat_Obj, npcs = 10)
  pca_map <- Embeddings(subset_Seurat_Obj, reduction = "pca")[rownames(subset_Seurat_Obj@meta.data),1:10]
  
  ### get slingshot object
  slingshot_obj <- slingshot(pca_map,
                             clusterLabels = subset_Seurat_Obj@meta.data$Development, 
                             reducedDim = "PCA")
  
  ### get colors for the clustering result
  cell_colors_clust <- cell_pal(unique(subset_Seurat_Obj@meta.data$Development), hue_pal())
  
  ### Trajectory inference
  png(paste0(outputDir2, "Trajectory_Inference_Without_Adult_PCA.png"), width = 2500, height = 1500, res = 200)
  plot(reducedDim(slingshot_obj),
       main=paste(type, "Trajectory Inference Without Adult"),
       col = cell_colors_clust[subset_Seurat_Obj@meta.data$Development],
       pch = 19, cex = 1)
  lines(slingshot_obj, lwd = 2, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         pch = 19)
  dev.off()
  
  ### Trajectory inference on multi dimentional PCA
  png(paste0(outputDir2, "Trajectory_Inference_Without_Adult_Multi-PCA.png"), width = 2500, height = 1500, res = 200)
  pairs(slingshot_obj, type="lineages", col = apply(slingshot_obj@clusterLabels, 1, function(x) cell_colors_clust[names(x)[which(x == 1)]]),
        show.constraints = TRUE, constraints.col = cell_colors_clust, cex = 0.8,
        horInd = 1:5, verInd = 1:5, main = paste0(type, "_Trajectory_Inference_Without_Adult"))
  par(xpd = TRUE)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         pch = 19, title = "Time")
  dev.off()
  
  ### DE & pathway analyses
  #
  # PC1
  multiple_analyses_in_one(Seurat_Object = subset_Seurat_Obj,
                           target_ident = NULL,
                           target_col = "Development",
                           target_col_factor_level = unique(subset_Seurat_Obj@meta.data$Development),
                           species = "mouse",
                           PC = "PC_1",
                           PC_Val = 0,
                           important_thresh = 0.1,
                           garbage_thresh = 1e-04,
                           result_dir = paste0(outputDir2, "PC1_0/"))
  
  
  ### new output directory
  type <- "LTHSC"
  outputDir2 <- paste0(outputDir, "ALL/", type, "/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### split the Seurat obj based on HSPC info
  subset_Seurat_Obj <- subset(Combined_Seurat_Obj, idents=type)
  
  ### order the meta data by developmental time
  subset_Seurat_Obj@meta.data <- subset_Seurat_Obj@meta.data[order(factor(subset_Seurat_Obj@meta.data$Development,
                                                                          levels = time_points)),]
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  subset_Seurat_Obj@assays$RNA@counts <- subset_Seurat_Obj@assays$RNA@counts[,rownames(subset_Seurat_Obj@meta.data)]
  
  ### run PCA
  subset_Seurat_Obj <- RunPCA(subset_Seurat_Obj, npcs = 10)
  pca_map <- Embeddings(subset_Seurat_Obj, reduction = "pca")[rownames(subset_Seurat_Obj@meta.data),1:10]
  
  ### get slingshot object
  slingshot_obj <- slingshot(pca_map,
                             clusterLabels = subset_Seurat_Obj@meta.data$Development, 
                             reducedDim = "PCA")
  
  ### get colors for the clustering result
  cell_colors_clust <- cell_pal(unique(subset_Seurat_Obj@meta.data$Development), hue_pal())
  
  ### Trajectory inference
  png(paste0(outputDir2, "Trajectory_Inference_Without_Adult_PCA.png"), width = 2500, height = 1500, res = 200)
  plot(reducedDim(slingshot_obj),
       main=paste(type, "Trajectory Inference Without Adult"),
       col = cell_colors_clust[subset_Seurat_Obj@meta.data$Development],
       pch = 19, cex = 1)
  lines(slingshot_obj, lwd = 2, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         pch = 19)
  dev.off()
  
  ### Trajectory inference on multi dimentional PCA
  png(paste0(outputDir2, "Trajectory_Inference_Without_Adult_Multi-PCA.png"), width = 2500, height = 1500, res = 200)
  pairs(slingshot_obj, type="lineages", col = apply(slingshot_obj@clusterLabels, 1, function(x) cell_colors_clust[names(x)[which(x == 1)]]),
        show.constraints = TRUE, constraints.col = cell_colors_clust, cex = 0.8,
        horInd = 1:5, verInd = 1:5, main = paste0(type, "_Trajectory_Inference_Without_Adult"))
  par(xpd = TRUE)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         pch = 19, title = "Time")
  dev.off()
  
  ### DE & pathway analyses
  #
  # PC1
  multiple_analyses_in_one(Seurat_Object = subset_Seurat_Obj,
                           target_ident = NULL,
                           target_col = "Development",
                           target_col_factor_level = unique(subset_Seurat_Obj@meta.data$Development),
                           species = "mouse",
                           PC = "PC_1",
                           PC_Val = -10,
                           important_thresh = 0.1,
                           garbage_thresh = 1e-04,
                           result_dir = paste0(outputDir2, "PC1_-10/"))
  
  
  #
  ### Analysis #1 & #2
  #
  
  ### remove ADULT cells
  subset_Seurat_Obj <- subset(Combined_Seurat_Obj, idents = unique_idents[-which(unique_idents == "ADULT")]) 
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  subset_Seurat_Obj@assays$RNA@counts <- subset_Seurat_Obj@assays$RNA@counts[,rownames(subset_Seurat_Obj@meta.data)]
  
  ### run PCA
  subset_Seurat_Obj <- FindVariableFeatures(subset_Seurat_Obj)
  subset_Seurat_Obj <- ScaleData(subset_Seurat_Obj)
  subset_Seurat_Obj <- RunPCA(subset_Seurat_Obj, npcs = 15)
  
  ### draw a PCA
  DimPlot(subset_Seurat_Obj, reduction = "pca", group.by = "Development", pt.size = 1.5) +
    labs(title = paste0("PCA_Subset"))
  
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
  
  ### Trajectory inference on multi dimentional PCA
  png(paste0(outputDir2, "Trajectory_Inference_Without_Adult_Multi-PCA.png"), width = 2500, height = 1500, res = 200)
  pairs(slingshot_obj, type="lineages", col = apply(slingshot_obj@clusterLabels, 1, function(x) cell_colors_clust[names(x)[which(x == 1)]]),
        show.constraints = TRUE, constraints.col = cell_colors_clust, cex = 0.8,
        horInd = 1:5, verInd = 1:5, main = paste0(type, "_Trajectory_Inference_Without_Adult"))
  par(xpd = TRUE)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         pch = 19, title = "Time")
  dev.off()
  
  ### DE & pathway analyses
  #
  # PC1
  multiple_analyses_in_one(Seurat_Object = subset_Seurat_Obj2,
                           target_ident = NULL,
                           target_col = "Development",
                           target_col_factor_level = unique(subset_Seurat_Obj2@meta.data$Development),
                           species = "mouse",
                           PC = "PC_1",
                           PC_Val = -10,
                           important_thresh = 0.1,
                           garbage_thresh = 1e-04,
                           result_dir = paste0(outputDir2, "PC1_-10/"))
  # PC2
  multiple_analyses_in_one(Seurat_Object = subset_Seurat_Obj2,
                           target_ident = NULL,
                           target_col = "Development",
                           target_col_factor_level = unique(subset_Seurat_Obj2@meta.data$Development),
                           species = "mouse",
                           PC = "PC_2",
                           PC_Val = 0,
                           important_thresh = 0.1,
                           garbage_thresh = 1e-04,
                           result_dir = paste0(outputDir2, "PC2_0/"))
  
  
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
  
  ### Trajectory inference on multi dimentional PCA
  png(paste0(outputDir2, "Trajectory_Inference_Without_Adult_Multi-PCA.png"), width = 2500, height = 1500, res = 200)
  pairs(slingshot_obj, type="lineages", col = apply(slingshot_obj@clusterLabels, 1, function(x) cell_colors_clust[names(x)[which(x == 1)]]),
        show.constraints = TRUE, constraints.col = cell_colors_clust, cex = 0.8,
        horInd = 1:5, verInd = 1:5, main = paste0(type, "_Trajectory_Inference_Without_Adult"))
  par(xpd = TRUE)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         pch = 19, title = "Time")
  dev.off()
  
  ### DE & pathway analyses
  #
  # PC1
  multiple_analyses_in_one(Seurat_Object = subset_Seurat_Obj2,
                           target_ident = NULL,
                           target_col = "Development",
                           target_col_factor_level = unique(subset_Seurat_Obj2@meta.data$Development),
                           species = "mouse",
                           PC = "PC_1",
                           PC_Val = -30,
                           important_thresh = 0.1,
                           garbage_thresh = 1e-04,
                           result_dir = paste0(outputDir2, "PC1_-30/"))
  # PC2
  multiple_analyses_in_one(Seurat_Object = subset_Seurat_Obj2,
                           target_ident = NULL,
                           target_col = "Development",
                           target_col_factor_level = unique(subset_Seurat_Obj2@meta.data$Development),
                           species = "mouse",
                           PC = "PC_2",
                           PC_Val = 10,
                           important_thresh = 0.1,
                           garbage_thresh = 1e-04,
                           result_dir = paste0(outputDir2, "PC2_10/"))
  # PC3
  multiple_analyses_in_one(Seurat_Object = subset_Seurat_Obj2,
                           target_ident = NULL,
                           target_col = "Development",
                           target_col_factor_level = unique(subset_Seurat_Obj2@meta.data$Development),
                           species = "mouse",
                           PC = "PC_3",
                           PC_Val = 10,
                           important_thresh = 0.1,
                           garbage_thresh = 1e-04,
                           result_dir = paste0(outputDir2, "PC3_10/"))
  
  
  #
  ### Analysis #3, #4, #5, & #6
  #
  
  ### updated stroma RDS file
  Updated_Seurat_Obj <- readRDS(file = Robj2_path)
  
  ### create new column for the analysis
  Updated_Seurat_Obj@meta.data$Dev_Anno <- paste0(Updated_Seurat_Obj@meta.data$Development, "_",
                                                  Updated_Seurat_Obj@meta.data$Annotation)
  
  ### determine necessary variables
  time_points <- c("E16", "E18", "P0", "ADULT")
  HSPC_populations <- c("LTHSC", "STHSC", "MPP2", "MPP3", "MPP4")
  
  ### set the ident of the object with the HSPC type
  Updated_Seurat_Obj <- SetIdent(object = Updated_Seurat_Obj,
                                 cells = rownames(Updated_Seurat_Obj@meta.data),
                                 value = Updated_Seurat_Obj@meta.data$Tissue)
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Updated_Seurat_Obj)), rownames(Updated_Seurat_Obj@meta.data)))
  
  
  ### new output directory
  type <- "StromaE16"
  outputDir2 <- paste0(outputDir, type, "/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### split the Seurat obj based on HSPC info
  subset_Seurat_Obj <- subset(Updated_Seurat_Obj, idents=type)
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  subset_Seurat_Obj@assays$RNA@counts <- subset_Seurat_Obj@assays$RNA@counts[,rownames(subset_Seurat_Obj@meta.data)]
  
  ### run PCA
  subset_Seurat_Obj <- RunPCA(subset_Seurat_Obj, npcs = 5)
  pca_map <- Embeddings(subset_Seurat_Obj, reduction = "pca")[rownames(subset_Seurat_Obj@meta.data),]
  
  ### get slingshot object
  slingshot_obj <- slingshot(pca_map,
                             clusterLabels = subset_Seurat_Obj@meta.data$Dev_Anno, 
                             reducedDim = "PCA")
  
  ### get colors for the clustering result
  cell_colors_clust <- cell_pal(unique(subset_Seurat_Obj@meta.data$Dev_Anno), hue_pal())
  
  ### Trajectory inference
  png(paste0(outputDir2, type, "_Trajectory_Inference_PCA.png"), width = 2500, height = 1500, res = 200)
  plot(reducedDim(slingshot_obj),
       main=paste(type, "Trajectory Inference (PCA)"),
       col = cell_colors_clust[as.character(subset_Seurat_Obj@meta.data$Dev_Anno)],
       pch = 19, cex = 1)
  lines(slingshot_obj, lwd = 2, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  dev.off()
  
  ### Trajectory inference on multi dimentional PCA
  png(paste0(outputDir2, type, "_Trajectory_Inference_Multi-PCA.png"), width = 2500, height = 1500, res = 200)
  pairs(slingshot_obj, type="lineages", col = apply(slingshot_obj@clusterLabels, 1, function(x) cell_colors_clust[names(x)[which(x == 1)]]),
        show.constraints = TRUE, constraints.col = cell_colors_clust, cex = 0.8,
        horInd = 1:5, verInd = 1:5, main = paste0(type, "_Trajectory_Inference"))
  par(xpd = TRUE)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         title = "Clusters",  pch = 19)
  dev.off()
  
  ### new output directory
  type <- "StromaE18"
  outputDir2 <- paste0(outputDir, type, "/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### split the Seurat obj based on HSPC info
  subset_Seurat_Obj <- subset(Updated_Seurat_Obj, idents=type)
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  subset_Seurat_Obj@assays$RNA@counts <- subset_Seurat_Obj@assays$RNA@counts[,rownames(subset_Seurat_Obj@meta.data)]
  
  ### run PCA
  subset_Seurat_Obj <- RunPCA(subset_Seurat_Obj, npcs = 5)
  pca_map <- Embeddings(subset_Seurat_Obj, reduction = "pca")[rownames(subset_Seurat_Obj@meta.data),]
  
  ### get slingshot object
  slingshot_obj <- slingshot(pca_map,
                             clusterLabels = subset_Seurat_Obj@meta.data$Dev_Anno, 
                             reducedDim = "PCA")
  
  ### get colors for the clustering result
  cell_colors_clust <- cell_pal(unique(subset_Seurat_Obj@meta.data$Dev_Anno), hue_pal())
  
  ### Trajectory inference
  png(paste0(outputDir2, type, "_Trajectory_Inference_PCA.png"), width = 2500, height = 1500, res = 200)
  plot(reducedDim(slingshot_obj),
       main=paste(type, "Trajectory Inference (PCA)"),
       col = cell_colors_clust[as.character(subset_Seurat_Obj@meta.data$Dev_Anno)],
       pch = 19, cex = 1)
  lines(slingshot_obj, lwd = 2, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  dev.off()
  
  ### Trajectory inference on multi dimentional PCA
  png(paste0(outputDir2, type, "_Trajectory_Inference_Multi-PCA.png"), width = 2500, height = 1500, res = 200)
  pairs(slingshot_obj, type="lineages", col = apply(slingshot_obj@clusterLabels, 1, function(x) cell_colors_clust[names(x)[which(x == 1)]]),
        show.constraints = TRUE, constraints.col = cell_colors_clust, cex = 0.8,
        horInd = 1:5, verInd = 1:5, main = paste0(type, "_Trajectory_Inference"))
  par(xpd = TRUE)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         title = "Clusters",  pch = 19)
  dev.off()
  
  ### new output directory
  type <- "StromaP0"
  outputDir2 <- paste0(outputDir, type, "/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### split the Seurat obj based on HSPC info
  subset_Seurat_Obj <- subset(Updated_Seurat_Obj, idents=type)
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  subset_Seurat_Obj@assays$RNA@counts <- subset_Seurat_Obj@assays$RNA@counts[,rownames(subset_Seurat_Obj@meta.data)]
  
  ### run PCA
  subset_Seurat_Obj <- RunPCA(subset_Seurat_Obj, npcs = 5)
  pca_map <- Embeddings(subset_Seurat_Obj, reduction = "pca")[rownames(subset_Seurat_Obj@meta.data),]
  
  ### get slingshot object
  slingshot_obj <- slingshot(pca_map,
                             clusterLabels = subset_Seurat_Obj@meta.data$Dev_Anno, 
                             reducedDim = "PCA")
  
  ### get colors for the clustering result
  cell_colors_clust <- cell_pal(unique(subset_Seurat_Obj@meta.data$Dev_Anno), hue_pal())
  
  ### Trajectory inference
  png(paste0(outputDir2, type, "_Trajectory_Inference_PCA.png"), width = 2500, height = 1500, res = 200)
  plot(reducedDim(slingshot_obj),
       main=paste(type, "Trajectory Inference (PCA)"),
       col = cell_colors_clust[as.character(subset_Seurat_Obj@meta.data$Dev_Anno)],
       pch = 19, cex = 1)
  lines(slingshot_obj, lwd = 2, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  dev.off()
  
  ### Trajectory inference on multi dimentional PCA
  png(paste0(outputDir2, type, "_Trajectory_Inference_Multi-PCA.png"), width = 2500, height = 1500, res = 200)
  pairs(slingshot_obj, type="lineages", col = apply(slingshot_obj@clusterLabels, 1, function(x) cell_colors_clust[names(x)[which(x == 1)]]),
        show.constraints = TRUE, constraints.col = cell_colors_clust, cex = 0.8,
        horInd = 1:5, verInd = 1:5, main = paste0(type, "_Trajectory_Inference"))
  par(xpd = TRUE)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         title = "Clusters",  pch = 19)
  dev.off()
  
  
  ### Analysis #5
  ### Pseudotime Analysis by combining all updated stroma objects,
  ### but for the adult library, only visualize the "CARs" population
  ### (visualize all other annotated populations from other libraries)
  
  ### new output directory
  type <- "Stroma_Adult_CARs"
  outputDir2 <- paste0(outputDir, type, "/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### set the ident of the object with the cell type to extract Stromas only
  Updated_Seurat_Obj <- SetIdent(object = Updated_Seurat_Obj,
                                 cells = rownames(Updated_Seurat_Obj@meta.data),
                                 value = Updated_Seurat_Obj@meta.data$Cell_Type)
  
  ### split the Seurat obj based on the cell type
  subset_Seurat_Obj <- subset(Updated_Seurat_Obj, idents="Stroma")
  
  ### create new column for the analysis
  subset_Seurat_Obj@meta.data$Dev_Anno <- paste0(subset_Seurat_Obj@meta.data$Development, "_",
                                                 subset_Seurat_Obj@meta.data$Annotation)
  
  ### set the ident of the object with the Trent's annotation
  subset_Seurat_Obj <- SetIdent(object = subset_Seurat_Obj,
                                cells = rownames(subset_Seurat_Obj@meta.data),
                                value = subset_Seurat_Obj@meta.data$Dev_Anno)
  
  ### split the Seurat obj based on the Trent's annotation
  subset_Seurat_Obj <- subset(subset_Seurat_Obj, idents=c("ADULT_CARs",
                                                          unique(subset_Seurat_Obj@meta.data$Dev_Anno)[!grepl("ADULT", unique(subset_Seurat_Obj@meta.data$Dev_Anno))]))
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  subset_Seurat_Obj@assays$RNA@counts <- subset_Seurat_Obj@assays$RNA@counts[,rownames(subset_Seurat_Obj@meta.data)]
  
  ### run UMAP
  subset_Seurat_Obj <- RunUMAP(subset_Seurat_Obj, dims = 1:5)
  umap_map <- Embeddings(subset_Seurat_Obj, reduction = "umap")[rownames(subset_Seurat_Obj@meta.data),]
  
  ### get slingshot object
  slingshot_obj <- slingshot(umap_map,
                             clusterLabels = subset_Seurat_Obj@meta.data$Development, 
                             reducedDim = "UMAP")
  
  ### get colors for the clustering result
  cell_colors_clust <- cell_pal(time_points, hue_pal())
  
  ### Trajectory inference
  png(paste0(outputDir2, type, "_Trajectory_Inference_UMAP.png"), width = 2500, height = 1500, res = 200)
  plot(reducedDim(slingshot_obj),
       main=paste(type, "Trajectory Inference (UMAP)"),
       col = cell_colors_clust[as.character(subset_Seurat_Obj@meta.data$Development)],
       pch = 19, cex = 1)
  lines(slingshot_obj, lwd = 2, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         title = "Clusters",  pch = 19)
  dev.off()
  
  ### run PCA
  subset_Seurat_Obj <- RunPCA(subset_Seurat_Obj, npcs = 5)
  pca_map <- Embeddings(subset_Seurat_Obj, reduction = "pca")[rownames(subset_Seurat_Obj@meta.data),]
  
  ### get slingshot object
  slingshot_obj <- slingshot(pca_map,
                             clusterLabels = subset_Seurat_Obj@meta.data$Development, 
                             reducedDim = "PCA")
  
  ### get colors for the clustering result
  cell_colors_clust <- cell_pal(time_points, hue_pal())
  
  ### Trajectory inference
  png(paste0(outputDir2, type, "_Trajectory_Inference_PCA.png"), width = 2500, height = 1500, res = 200)
  plot(reducedDim(slingshot_obj),
       main=paste(type, "Trajectory Inference (PCA)"),
       col = cell_colors_clust[as.character(subset_Seurat_Obj@meta.data$Development)],
       pch = 19, cex = 1)
  lines(slingshot_obj, lwd = 2, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         title = "Clusters",  pch = 19)
  dev.off()
  
  ### Trajectory inference on multi dimentional PCA
  png(paste0(outputDir2, type, "_Trajectory_Inference_Multi-PCA.png"), width = 2500, height = 1500, res = 200)
  pairs(slingshot_obj, type="lineages", col = apply(slingshot_obj@clusterLabels, 1, function(x) cell_colors_clust[names(x)[which(x == 1)]]),
        show.constraints = TRUE, constraints.col = cell_colors_clust, cex = 0.8,
        horInd = 1:5, verInd = 1:5, main = paste0(type, "_Trajectory_Inference"))
  par(xpd = TRUE)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         title = "Clusters",  pch = 19)
  dev.off()
  
  
  ### Analysis #6
  ### Perform with all of the updated Heme objects.
  ### only visualize the "T Cell" and "Mast" clusters in E16.5,
  ### along with all other clusters from the other libraries
  
  ### new output directory
  type <- "Heme_E16_TCell_Mast"
  outputDir2 <- paste0(outputDir, type, "/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### set the ident of the object with the cell type to extract Stromas only
  Updated_Seurat_Obj <- SetIdent(object = Updated_Seurat_Obj,
                                 cells = rownames(Updated_Seurat_Obj@meta.data),
                                 value = Updated_Seurat_Obj@meta.data$Cell_Type)
  
  ### split the Seurat obj based on the cell type
  subset_Seurat_Obj <- subset(Updated_Seurat_Obj, idents="Heme")
  
  ### create new column for the analysis
  subset_Seurat_Obj@meta.data$Dev_Anno <- paste0(subset_Seurat_Obj@meta.data$Development, "_",
                                                 subset_Seurat_Obj@meta.data$Annotation)
  
  ### set the ident of the object with the Trent's annotation
  subset_Seurat_Obj <- SetIdent(object = subset_Seurat_Obj,
                                cells = rownames(subset_Seurat_Obj@meta.data),
                                value = subset_Seurat_Obj@meta.data$Dev_Anno)
  
  ### split the Seurat obj based on the Trent's annotation
  subset_Seurat_Obj <- subset(subset_Seurat_Obj, idents=c("E16_Mono/Mast", "E16_Mast1", "E16_Mast2", "E16_Mast3",
                                                          "E16_TCell1", "E16_TCell2", "E16_TCell3", "E16_TCell4", "E16_TCell5", "E16_TCell6",
                                                          unique(subset_Seurat_Obj@meta.data$Dev_Anno)[!grepl("E16", unique(subset_Seurat_Obj@meta.data$Dev_Anno))]))
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  subset_Seurat_Obj@assays$RNA@counts <- subset_Seurat_Obj@assays$RNA@counts[,rownames(subset_Seurat_Obj@meta.data)]
  
  ### run UMAP
  subset_Seurat_Obj <- RunUMAP(subset_Seurat_Obj, dims = 1:5)
  umap_map <- Embeddings(subset_Seurat_Obj, reduction = "umap")[rownames(subset_Seurat_Obj@meta.data),]
  
  ### get slingshot object
  slingshot_obj <- slingshot(umap_map,
                             clusterLabels = subset_Seurat_Obj@meta.data$Dev_Anno, 
                             reducedDim = "UMAP")
  
  ### get colors for the clustering result
  cell_colors_clust <- cell_pal(unique(subset_Seurat_Obj@meta.data$Dev_Anno), hue_pal())
  
  ### Trajectory inference
  png(paste0(outputDir2, type, "_Trajectory_Inference_UMAP.png"), width = 2500, height = 1500, res = 200)
  plot(reducedDim(slingshot_obj),
       main=paste(type, "Trajectory Inference (UMAP)"),
       col = cell_colors_clust[as.character(subset_Seurat_Obj@meta.data$Dev_Anno)],
       pch = 19, cex = 1)
  lines(slingshot_obj, lwd = 2, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  dev.off()
  
  ### run PCA
  subset_Seurat_Obj <- RunPCA(subset_Seurat_Obj, npcs = 5)
  pca_map <- Embeddings(subset_Seurat_Obj, reduction = "pca")[rownames(subset_Seurat_Obj@meta.data),]
  
  ### get slingshot object
  slingshot_obj <- slingshot(pca_map,
                             clusterLabels = subset_Seurat_Obj@meta.data$Dev_Anno, 
                             reducedDim = "PCA")
  
  ### get colors for the clustering result
  cell_colors_clust <- cell_pal(unique(subset_Seurat_Obj@meta.data$Dev_Anno), hue_pal())
  
  ### Trajectory inference
  png(paste0(outputDir2, type, "_Trajectory_Inference_PCA.png"), width = 2500, height = 1500, res = 200)
  plot(reducedDim(slingshot_obj),
       main=paste(type, "Trajectory Inference (PCA)"),
       col = cell_colors_clust[as.character(subset_Seurat_Obj@meta.data$Dev_Anno)],
       pch = 19, cex = 1)
  lines(slingshot_obj, lwd = 2, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  dev.off()
  
  ### Trajectory inference on multi dimentional PCA
  png(paste0(outputDir2, type, "_Trajectory_Inference_Multi-PCA.png"), width = 2500, height = 1500, res = 200)
  pairs(slingshot_obj, type="lineages", col = apply(slingshot_obj@clusterLabels, 1, function(x) cell_colors_clust[names(x)[which(x == 1)]]),
        show.constraints = TRUE, constraints.col = cell_colors_clust, cex = 0.8,
        horInd = 1:5, verInd = 1:5, main = paste0(type, "_Trajectory_Inference"))
  par(xpd = TRUE)
  dev.off()
  
  
  #
  ### RNA Magnet
  #
  
  ### merge specific Heme and Stroma and run RNA-Magnet
  merge_heme_stroma_and_run_rnamagnet <- function(Seurat_Object,
                                                  target_col,
                                                  comp1,
                                                  comp2,
                                                  time_point,
                                                  dim_method=c("UMAP", "PCA"),
                                                  result_dir="./") {
    
    ### set output directory
    result_dir <- paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                         "_RNAMagnet_Result/", time_point, "/")
    dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
    
    ### set group info to the metadata
    Seurat_Object@meta.data$Group <- paste0(Seurat_Object@meta.data[,target_col], "_", Seurat_Object@meta.data$Development)
    
    ### all the comps
    comps <- union(paste0(comp1, "_", time_point), paste0(comp2, "_", time_point))
    
    ### set the ident of the object with the specified info
    Seurat_Object <- SetIdent(object = Seurat_Object,
                              cells = rownames(Seurat_Object@meta.data),
                              value = Seurat_Object@meta.data$Group)
    
    ### only keep the specified cells
    Seurat_Object <- subset(Seurat_Object, idents=comps)
    
    ### check whether the orders are the same
    print(identical(names(Idents(object = Seurat_Object)), rownames(Seurat_Object@meta.data)))
    
    ### rownames in the meta.data should be in the same order as colnames in the counts
    Seurat_Object@meta.data <- Seurat_Object@meta.data[colnames(Seurat_Object@assays$RNA@counts),]
    
    ### test
    # Seurat_Object <- subset(Seurat_Object, cells = rownames(Seurat_Object@meta.data)[which(Seurat_Object@reductions$pca@cell.embeddings[,"PC_2"] < 0)])
    
    ### preprocessing
    Seurat_Object <- FindVariableFeatures(Seurat_Object)
    Seurat_Object <- ScaleData(Seurat_Object)
    
    ### run PCA/UMAP
    if(dim_method == "PCA") {
      Seurat_Object <- RunPCA(Seurat_Object, npcs = 15)
      dim_map <- Embeddings(Seurat_Object, reduction = "pca")[rownames(Seurat_Object@meta.data),]
    } else if(dim_method == "UMAP") {
      Seurat_Object <- RunUMAP(Seurat_Object, dims = 1:15)
      dim_map <- Embeddings(Seurat_Object, reduction = "umap")[rownames(Seurat_Object@meta.data),]
    } else {
      stop("ERROR: dim_method not PCA nor UMAP.")
    }
    
    ### run RNAMagnet with anchors
    result <- RNAMagnetAnchors(Seurat_Object,
                               anchors = unique(Seurat_Object@meta.data$Group))
    
    ### write the result as an Excel file
    write.xlsx2(data.frame(Cell=rownames(result), result,
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                              "_RNAMagnet_Result_", time_point, ".xlsx"),
                sheetName = paste0("RNAMagnet_Result"),
                row.names = FALSE)
    
    ### add RNAMagnet info to the seurat object
    Seurat_Object@meta.data$direction <- as.character(result[rownames(Seurat_Object@meta.data),"direction"])
    Seurat_Object@meta.data$adhesiveness <- as.numeric(result[rownames(Seurat_Object@meta.data),"adhesiveness"])
    Seurat_Object@meta.data$specificity <- as.numeric(sapply(rownames(Seurat_Object@meta.data), function(x) {
      result[x, result[x,"direction"]]
    }))
    
    ### check the order
    print(identical(rownames(Seurat_Object@meta.data), rownames(dim_map)))
    
    ### make a data frame for ggplot
    plot_df <- data.frame(X=dim_map[rownames(Seurat_Object@meta.data),1],
                          Y=dim_map[rownames(Seurat_Object@meta.data),2],
                          direction = Seurat_Object@meta.data$direction,
                          group_alpha = Seurat_Object@meta.data$adhesiveness,
                          cluster_color = Seurat_Object@meta.data$Group,
                          specificity = Seurat_Object@meta.data$specificity,
                          stringsAsFactors = FALSE, check.names = FALSE)
    
    ### get colors for the clustering result
    cell_colors_clust <- cell_pal(unique(Seurat_Object@meta.data$Group), hue_pal())
    
    ### scatter plot
    p <- list()
    
    ### draw a scatter plot with the adhesiveness info
    p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Cluster") +
      ggtitle(paste0(dim_method, "with Cell Type")) +
      scale_color_brewer(palette="Dark2") +
      theme_classic(base_size = 16)
    
    p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="direction", alpha="group_alpha"), size=2) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Direction", alpha="Adhesiveness") +
      ggtitle(paste(dim_method, "with Direction & Adhesiveness")) +
      theme_classic(base_size = 16) +
      scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$direction)[order(unique(plot_df$direction))])],
                         labels = names(cell_colors_clust[as.character(unique(plot_df$direction)[order(unique(plot_df$direction))])]))
    
    ### save the plots
    g <- arrangeGrob(grobs = p,
                     nrow = 2,
                     ncol = 1,
                     top = "")
    ggsave(file = paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                         "_RNAMagnet_Result_AD_", time_point, ".png"), g, width = 20, height = 12, dpi = 300)
    
    ### draw a beeswarm plot with the adhesiveness info
    ggplot(plot_df, aes_string(x="cluster_color", y="group_alpha")) +
      theme_classic(base_size = 16) +
      geom_boxplot() +
      geom_beeswarm(aes_string(color="direction"), na.rm = TRUE) +
      stat_compare_means() +
      labs(x = "", y = "Adhesiveness") +
      theme(legend.position="right")
    ggsave(file = paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                         "_RNAMagnet_Beeswarm_AD_", time_point, ".png"), width = 20, height = 12, dpi = 300)
    
    ### draw a scatter plot with the specificity info
    p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Cluster") +
      ggtitle(paste(dim_method, "with Cell Type")) +
      scale_color_brewer(palette="Dark2") +
      theme_classic(base_size = 16)
    
    p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="direction", alpha="specificity"), size=2) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Direction", alpha="Specificity Score") +
      ggtitle(paste(dim_method, "with Direction & Specificity Score")) +
      theme_classic(base_size = 16) +
      scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$direction)[order(unique(plot_df$direction))])],
                         labels = names(cell_colors_clust[as.character(unique(plot_df$direction)[order(unique(plot_df$direction))])]))
    
    ### save the plots
    g <- arrangeGrob(grobs = p,
                     nrow = 2,
                     ncol = 1,
                     top = "")
    ggsave(file = paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                         "_RNAMagnet_Result_SP_", time_point, ".png"), g, width = 20, height = 12, dpi = 300)
    
    ### draw a beeswarm plot with the adhesiveness info
    ggplot(plot_df, aes_string(x="cluster_color", y="specificity")) +
      theme_classic(base_size = 16) +
      geom_boxplot() +
      geom_beeswarm(aes_string(color="direction"), na.rm = TRUE) +
      stat_compare_means() +
      labs(x = "", y = "Specificity") +
      theme(legend.position="right")
    ggsave(file = paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                         "_RNAMagnet_Beeswarm_SP_", time_point, ".png"), width = 20, height = 12, dpi = 300)
    
    
    ### run RNAMagnet signaling
    result2 <- RNAMagnetSignaling(Seurat_Object)
    
    ### draw signaling network
    png(paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
               "_Signaling_RNAMagnet_Network_", time_point, ".png"), width = 2400, height = 1200, res = 120)
    set.seed(1234)
    PlotSignalingNetwork(result2, threshold = 0.01,
                         colors = cell_colors_clust,
                         useLabels = TRUE,
                         main = paste0("Signaling_RNAMagnet_", time_point))
    legend("left",
           title = "Clusters",
           legend = names(cell_colors_clust),
           fill = cell_colors_clust,
           border = NA,
           bty = "o",
           x.intersp = 1,
           y.intersp = 0.8,
           cex = 1)
    dev.off()
    
    ### get all the interaction list
    interaction_list <- NULL
    for(clust1 in names(cell_colors_clust)) {
      for(clust2 in names(cell_colors_clust)) {
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
                file = paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                              "_Signaling_RNAMagnet_Interaction_List_", time_point, ".xlsx"),
                sheetName = paste0("Signaling_RNAMagnet_Interaction_List"),
                row.names = FALSE)
    
  }
  
  ### result directory
  outputDir2 <- paste0(outputDir, "RNA_Magnet/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### see python environment since RNAMagnet uses the python module 'magic'
  # conda_create("r-reticulate")
  # conda_install("r-reticulate", "python-magic")
  use_condaenv("r-reticulate")
  py_config()
  py_module_available("magic")
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Updated_Seurat_Obj)), rownames(Updated_Seurat_Obj@meta.data)))
  
  ### run UMAP
  Updated_Seurat_Obj <- FindVariableFeatures(Updated_Seurat_Obj)
  Updated_Seurat_Obj <- ScaleData(Updated_Seurat_Obj)
  Updated_Seurat_Obj <- RunUMAP(Updated_Seurat_Obj, dims = 1:5)
  
  ### draw a UMAP with cell type
  umap_plot <- DimPlot(Updated_Seurat_Obj, reduction = "umap", group.by = "Cell_Type", pt.size = 1.5) +
    labs(title = paste0("UMAP_Combined_Tissue"))
  umap_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.3
  plot(umap_plot)
  
  ### set the ident of the object with the group info
  Updated_Seurat_Obj <- SetIdent(object = Updated_Seurat_Obj,
                                   cells = rownames(Updated_Seurat_Obj@meta.data),
                                   value = Updated_Seurat_Obj@meta.data$Development)
  
  ### split the Seurat obj based on the given info
  Combined_Adult_Seurat_Obj <- subset(Updated_Seurat_Obj, idents="ADULT")
  
  ### draw a UMAP with cell type
  umap_plot <- DimPlot(Combined_Adult_Seurat_Obj, reduction = "umap", group.by = "Cell_Type", pt.size = 1.5) +
    labs(title = paste0("UMAP_Combined_Adults_Only"))
  umap_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  plot(umap_plot)
  ggsave(file = paste0(outputDir2, "UMAP_Combined_Adult_Cells.png"), width = 15, height = 10, dpi = 300)
  
  ### draw a UMAP with HSPC info
  umap_plot <- DimPlot(Combined_Adult_Seurat_Obj, reduction = "umap", group.by = "HSPC", pt.size = 1.5) +
    labs(title = paste0("UMAP_Combined_Adults_Only"))
  umap_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  plot(umap_plot)
  ggsave(file = paste0(outputDir2, "UMAP_Combined_Adult_HSPC.png"), width = 15, height = 10, dpi = 300)
  
  ### draw a UMAP with Trent's annotation
  Combined_Adult_Seurat_Obj@meta.data$Annotation <- factor(Combined_Adult_Seurat_Obj@meta.data$Annotation)
  umap_plot <- DimPlot(Combined_Adult_Seurat_Obj, reduction = "umap", group.by = "Annotation", pt.size = 1.5) +
    labs(title = paste0("UMAP_Combined_Adult_With_the_Annotation"))
  umap_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  umap_plot <- LabelClusters(plot = umap_plot, id = "Annotation", col = "black")
  ggsave(file = paste0(outputDir2, "UMAP_Combined_Adult_With_the_Annotation.png"), width = 15, height = 10, dpi = 300)
  
  ### LTHSC vs Stroma - Adult
  merge_heme_stroma_and_run_rnamagnet(Seurat_Object = Combined_Adult_Seurat_Obj,
                                      target_col = "HSPC",
                                      comp1 = "LTHSC",
                                      comp2 = "Stroma",
                                      time_point = "Adult",
                                      dim_method = "UMAP",
                                      result_dir=outputDir2)
  
  ### (STHSC, MPP2, MPP3, MPP4) vs Stroma - Adult
  merge_heme_stroma_and_run_rnamagnet(Seurat_Object = Combined_Adult_Seurat_Obj,
                                      target_col = "HSPC",
                                      comp1 = c("STHSC", "MPP2", "MPP3", "MPP4"),
                                      comp2 = "Stroma",
                                      time_point = "Adult",
                                      dim_method = "UMAP",
                                      result_dir=outputDir2)
  
  
  `%||%` <- function(lhs, rhs) {
    if (!is.null(x = lhs)) {
      return(lhs)
    } else {
      return(rhs)
    }
  }
  
  #' Label clusters on a ggplot2-based scatter plot
  #'
  #' @param plot A ggplot2-based scatter plot
  #' @param id Name of variable used for coloring scatter plot
  #' @param clusters Vector of cluster ids to label
  #' @param labels Custom labels for the clusters
  #' @param split.by Split labels by some grouping label, useful when using
  #' \code{\link[ggplot2]{facet_wrap}} or \code{\link[ggplot2]{facet_grid}}
  #' @param repel Use \code{geom_text_repel} to create nicely-repelled labels
  #' @param geom Name of geom to get X/Y aesthetic names for
  #' @param box Use geom_label/geom_label_repel (includes a box around the text
  #' labels)
  #' @param position How to place the label if repel = FALSE. If "median", place
  #' the label at the median position. If "nearest" place the label at the
  #' position of the nearest data point to the median.
  #' @param ... Extra parameters to \code{\link[ggrepel]{geom_text_repel}}, such as \code{size}
  #'
  #' @return A ggplot2-based scatter plot with cluster labels
  #'
  #' @importFrom stats median
  #' @importFrom ggrepel geom_text_repel geom_label_repel
  #' @importFrom ggplot2 aes_string geom_text geom_label
  #' @importFrom RANN nn2
  #'
  #' @export
  #'
  #' @seealso \code{\link[ggrepel]{geom_text_repel}} \code{\link[ggplot2]{geom_text}}
  #'
  #' @examples
  #' plot <- DimPlot(object = pbmc_small)
  #' LabelClusters(plot = plot, id = 'ident')
  #'
  LabelClusters <- function(
    plot,
    id,
    clusters = NULL,
    labels = NULL,
    split.by = NULL,
    repel = TRUE,
    box = FALSE,
    geom = 'GeomPoint',
    position = "median",
    ...
  ) {
    xynames <- unlist(x = GetXYAesthetics(plot = plot, geom = geom), use.names = TRUE)
    if (!id %in% colnames(x = plot$data)) {
      stop("Cannot find variable ", id, " in plotting data")
    }
    if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
      warning("Cannot find splitting variable ", id, " in plotting data")
      split.by <- NULL
    }
    data <- plot$data[, c(xynames, id, split.by)]
    possible.clusters <- as.character(x = na.omit(object = unique(x = data[, id])))
    groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[, id])))
    if (any(!groups %in% possible.clusters)) {
      stop("The following clusters were not found: ", paste(groups[!groups %in% possible.clusters], collapse = ","))
    }
    pb <- ggplot_build(plot = plot)
    if (geom == 'GeomSpatial') {
      data[, xynames["y"]] = max(data[, xynames["y"]]) - data[, xynames["y"]] + min(data[, xynames["y"]])
      if (!pb$plot$plot_env$crop) {
        y.transform <- c(0, nrow(x = pb$plot$plot_env$image)) - pb$layout$panel_params[[1]]$y.range
        data[, xynames["y"]] <- data[, xynames["y"]] + sum(y.transform)
      }
    }
    data <- cbind(data, color = pb$data[[1]][[1]])
    labels.loc <- lapply(
      X = groups,
      FUN = function(group) {
        data.use <- data[data[, id] == group, , drop = FALSE]
        data.medians <- if (!is.null(x = split.by)) {
          do.call(
            what = 'rbind',
            args = lapply(
              X = unique(x = data.use[, split.by]),
              FUN = function(split) {
                medians <- apply(
                  X = data.use[data.use[, split.by] == split, xynames, drop = FALSE],
                  MARGIN = 2,
                  FUN = median,
                  na.rm = TRUE
                )
                medians <- as.data.frame(x = t(x = medians))
                medians[, split.by] <- split
                return(medians)
              }
            )
          )
        } else {
          as.data.frame(x = t(x = apply(
            X = data.use[, xynames, drop = FALSE],
            MARGIN = 2,
            FUN = median,
            na.rm = TRUE
          )))
        }
        data.medians[, id] <- group
        data.medians$color <- data.use$color[1]
        return(data.medians)
      }
    )
    if (position == "nearest") {
      labels.loc <- lapply(X = labels.loc, FUN = function(x) {
        group.data <- data[as.character(x = data[, id]) == as.character(x[3]), ]
        nearest.point <- nn2(data = group.data[, 1:2], query = as.matrix(x = x[c(1,2)]), k = 1)$nn.idx
        x[1:2] <- group.data[nearest.point, 1:2]
        return(x)
      })
    }
    labels.loc <- do.call(what = 'rbind', args = labels.loc)
    if(is.null(levels(data[, id]))) {
      labels.loc[, id] <- factor(x = labels.loc[, id], levels = unique(data[, id]))
    } else {
      labels.loc[, id] <- factor(x = labels.loc[, id], levels = levels(data[, id]))
    }
    labels <- labels %||% groups
    if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
      stop("Length of labels (", length(x = labels),  ") must be equal to the number of clusters being labeled (", length(x = labels.loc), ").")
    }
    names(x = labels) <- groups
    for (group in groups) {
      labels.loc[labels.loc[, id] == group, id] <- labels[group]
    }
    if (box) {
      geom.use <- ifelse(test = repel, yes = geom_label_repel, no = geom_label)
      plot <- plot + geom.use(
        data = labels.loc,
        mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id, fill = id),
        show.legend = FALSE,
        ...
      ) + scale_fill_manual(values = labels.loc$color[order(labels.loc[, id])])
    } else {
      geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
      plot <- plot + geom.use(
        data = labels.loc,
        mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
        show.legend = FALSE,
        ...
      )
    }
    return(plot)
  }
  
  
  ### DE between two directions in Stroma
  Seurat_Object = Combined_Adult_Seurat_Obj
  target_col = "HSPC"
  comp1 = "LTHSC"
  comp2 = "Stroma"
  time_point = "Adult"
  dim_method = "UMAP"
  result_dir=outputDir2
  
  ### set output directory
  result_dir <- paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                       "_RNAMagnet_Result/", time_point, "/")
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
  
  ### set group info to the metadata
  Seurat_Object@meta.data$Group <- paste0(Seurat_Object@meta.data[,target_col], "_", Seurat_Object@meta.data$Development)
  
  ### all the comps
  comps <- union(paste0(comp1, "_", time_point), paste0(comp2, "_", time_point))
  
  ### set the ident of the object with the specified info
  Seurat_Object <- SetIdent(object = Seurat_Object,
                            cells = rownames(Seurat_Object@meta.data),
                            value = Seurat_Object@meta.data$Group)
  
  ### only keep the specified cells
  Seurat_Object <- subset(Seurat_Object, idents=comps)
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Object)), rownames(Seurat_Object@meta.data)))
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Object@meta.data <- Seurat_Object@meta.data[colnames(Seurat_Object@assays$RNA@counts),]
  
  ### preprocessing
  Seurat_Object <- FindVariableFeatures(Seurat_Object)
  Seurat_Object <- ScaleData(Seurat_Object)
  
  ### run PCA/UMAP
  if(dim_method == "PCA") {
    Seurat_Object <- RunPCA(Seurat_Object, npcs = 15)
    dim_map <- Embeddings(Seurat_Object, reduction = "pca")[rownames(Seurat_Object@meta.data),]
  } else if(dim_method == "UMAP") {
    Seurat_Object <- RunUMAP(Seurat_Object, dims = 1:15)
    dim_map <- Embeddings(Seurat_Object, reduction = "umap")[rownames(Seurat_Object@meta.data),]
  } else {
    stop("ERROR: dim_method not PCA nor UMAP.")
  }
  
  ### draw a UMAP with Trent's annotation
  Seurat_Object@meta.data$Annotation <- factor(Seurat_Object@meta.data$Annotation)
  umap_plot <- DimPlot(Seurat_Object, reduction = "umap", group.by = "Annotation", pt.size = 1.5) +
    labs(title = paste0("UMAP_Combined_Adult_With_the_Annotation"))
  umap_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  umap_plot <- LabelClusters(plot = umap_plot, id = "Annotation", col = "black")
  ggsave(file = paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                       "_UMAP_with_the_annotation_", time_point, ".png"),
         width = 15, height = 10, dpi = 300)
  
  ### run RNAMagnet with anchors
  ### Warning: the Idents(Seurat_Object) should be along with the given 'anchors' input
  ### Another ERROR:
  ### The following code is fixed due to an ERROR: 'drop=FALSE'
  #compute mean gene expression level per population
  # out@anchors <- do.call(cbind, lapply(anchors, function(id) {
  #   apply(resolvedRawData[,seurat.ident == id,drop=FALSE],1,mean)
  # }));
  # colnames(out@anchors) <- anchors
  Seurat_Object <- SetIdent(object = Seurat_Object,
                            cells = rownames(Seurat_Object@meta.data),
                            value = Seurat_Object@meta.data$Group)
  result <- RNAMagnetAnchors(Seurat_Object,
                             anchors = unique(Seurat_Object@meta.data$Group))
  Seurat_Object <- SetIdent(object = Seurat_Object,
                            cells = rownames(Seurat_Object@meta.data),
                            value = Seurat_Object@meta.data$Annotation)
  result2 <- RNAMagnetAnchors(Seurat_Object,
                              anchors = unique(Seurat_Object@meta.data$Annotation))
  
  ### add RNAMagnet info to the seurat object
  Seurat_Object@meta.data$direction <- as.character(result[rownames(Seurat_Object@meta.data),"direction"])
  Seurat_Object@meta.data$direction2 <- as.character(result2[rownames(Seurat_Object@meta.data),"direction"])
  Seurat_Object@meta.data$adhesiveness <- as.numeric(result[rownames(Seurat_Object@meta.data),"adhesiveness"])
  Seurat_Object@meta.data$adhesiveness2 <- as.numeric(result2[rownames(Seurat_Object@meta.data),"adhesiveness"])
  Seurat_Object@meta.data$specificity <- as.numeric(sapply(rownames(Seurat_Object@meta.data), function(x) {
    result[x, result[x,"direction"]]
  }))
  Seurat_Object@meta.data$specificity2 <- as.numeric(sapply(rownames(Seurat_Object@meta.data), function(x) {
    result2[x, result2[x,"direction"]]
  }))
  
  ### check the order
  print(identical(rownames(Seurat_Object@meta.data), rownames(dim_map)))
  
  ### make a data frame for ggplot
  plot_df <- data.frame(X=dim_map[rownames(Seurat_Object@meta.data),1],
                        Y=dim_map[rownames(Seurat_Object@meta.data),2],
                        direction = Seurat_Object@meta.data$direction,
                        direction2 = Seurat_Object@meta.data$direction2,
                        adhesiveness = Seurat_Object@meta.data$adhesiveness,
                        adhesiveness2 = Seurat_Object@meta.data$adhesiveness2,
                        cluster_color = Seurat_Object@meta.data$Group,
                        specificity = Seurat_Object@meta.data$specificity,
                        specificity2 = Seurat_Object@meta.data$specificity2,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### get colors for the clustering result
  cell_colors_clust <- cell_pal(unique(Seurat_Object@meta.data$Annotation), hue_pal())
  
  ### scatter plot
  p <- list()
  
  ### draw a scatter plot with the adhesiveness info
  p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
    xlab("PC1") + ylab("PC2") +
    labs(col="Cluster") +
    ggtitle("UMAP with Cell Type") +
    scale_color_brewer(palette="Dark2") +
    theme_classic(base_size = 16)
  
  p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="direction", alpha="adhesiveness"), size=2) +
    xlab("PC1") + ylab("PC2") +
    labs(col="Direction", alpha="Adhesiveness") +
    ggtitle("UMAP with Direction & Adhesiveness") +
    scale_color_brewer(palette="Set1") +
    theme_classic(base_size = 16)
  
  plot_df$direction2 <- factor(plot_df$direction2, levels = unique(plot_df$direction2))
  p[[3]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="direction2", alpha="adhesiveness2"), size=2) +
    xlab("PC1") + ylab("PC2") +
    labs(col="Direction", alpha="Adhesiveness") +
    ggtitle("UMAP with Direction & Adhesiveness") +
    theme_classic(base_size = 16) +
    scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])],
                       labels = names(cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])]))
  ### the id column in the plot_df should be a factor
  p[[3]] <- LabelClusters(plot = p[[3]], id = "direction2", col = "black")
  
  ### save the plots
  g <- arrangeGrob(grobs = p,
                   nrow = 3,
                   ncol = 1,
                   top = "")
  ggsave(file = paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                       "_RNAMagnet_Result_AD_", time_point, "_New_Annotation.png"), g, width = 20, height = 20, dpi = 300)
  
  ### draw a scatter plot with the specificity info
  p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="cluster_color"), size=2, alpha=0.5) +
    xlab("PC1") + ylab("PC2") +
    labs(col="Cluster") +
    ggtitle("UMAP with Cell Type") +
    scale_color_brewer(palette="Dark2") +
    theme_classic(base_size = 16)
  
  p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="direction", alpha="specificity"), size=2) +
    xlab("PC1") + ylab("PC2") +
    labs(col="Direction", alpha="Specificity") +
    ggtitle("UMAP with Direction & Specificity") +
    scale_color_brewer(palette="Set1") +
    theme_classic(base_size = 16)
  
  plot_df$direction2 <- factor(plot_df$direction2, levels = unique(plot_df$direction2))
  p[[3]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
    geom_point(aes_string(col="direction2", alpha="specificity2"), size=2) +
    xlab("PC1") + ylab("PC2") +
    labs(col="Direction", alpha="Specificity") +
    ggtitle("UMAP with Direction & Specificity") +
    theme_classic(base_size = 16) +
    scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])],
                       labels = names(cell_colors_clust[as.character(unique(plot_df$direction2)[order(unique(plot_df$direction2))])]))
  ### the id column in the plot_df should be a factor
  p[[3]] <- LabelClusters(plot = p[[3]], id = "direction2", col = "black")
  
  ### save the plots
  g <- arrangeGrob(grobs = p,
                   nrow = 3,
                   ncol = 1,
                   top = "")
  ggsave(file = paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                       "_RNAMagnet_Result_SP_", time_point, "_New_Annotation.png"), g, width = 20, height = 20, dpi = 300)
  
  
  ### set the ident of the object with the Cell Type
  Seurat_Object <- SetIdent(object = Seurat_Object,
                            cells = rownames(Seurat_Object@meta.data),
                            value = Seurat_Object@meta.data$Cell_Type)
  
  ### only keep the specified cells
  Seurat_Object <- subset(Seurat_Object, idents="Stroma")
  
  ### set the ident of the object with the Cell Type
  Seurat_Object <- SetIdent(object = Seurat_Object,
                            cells = rownames(Seurat_Object@meta.data),
                            value = Seurat_Object@meta.data$direction)
  
  ### get DE genes between two groups
  de_result <- FindMarkers(object = Seurat_Object, test.use = "wilcox",
                           ident.1 = "LTHSC_Adult", ident.2 = "Stroma_Adult",
                           logfc.threshold = 0, min.pct = 0.1)
  
  ### order the DE reuslt
  de_result <- de_result[order(de_result$p_val_adj),]
  
  ### data frame
  de_result <- data.frame(Gene_Symbol=rownames(de_result),
                          de_result,
                          stringsAsFactors = FALSE, check.names = FALSE)
  
  ### write out the DE result
  write.xlsx2(de_result,
              file = paste0(result_dir,
                            paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                            "_RNAMagnet_DE_Result_", time_point, ".xlsx"),
              sheetName = "DE_Result",
              row.names = FALSE)
  
  ### pathway analysis
  pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Mm.eg.db,
                                                            de_result[which(de_result$p_val_adj < 1e-100),
                                                                      "Gene_Symbol"],
                                                            "ENTREZID", "SYMBOL"),
                                          org = "mouse", database = "GO",
                                          title = paste0(paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                                                         "_RNAMagnet_Pathway_Result"),
                                          displayNum = 30, imgPrint = TRUE,
                                          dir = paste0(result_dir))
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Mm.eg.db,
                                                              de_result[which(de_result$p_val_adj < 1e-100),
                                                                        "Gene_Symbol"],
                                                              "ENTREZID", "SYMBOL"),
                                            org = "mouse", database = "KEGG",
                                            title = paste0(paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                                                           "_RNAMagnet_Pathway_Result"),
                                            displayNum = 30, imgPrint = TRUE,
                                            dir = paste0(result_dir))
  if(!is.null(pathway_result_GO) && nrow(pathway_result_GO) > 0) {
    write.xlsx2(pathway_result_GO, file = paste0(result_dir,
                                                 paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                                                 "_RNAMagnet_GO_Result_", time_point, ".xlsx"),
                row.names = FALSE, sheetName = paste0("GO_Results"))
  }
  if(!is.null(pathway_result_KEGG) && nrow(pathway_result_KEGG) > 0) {
    write.xlsx2(pathway_result_KEGG, file = paste0(result_dir,
                                                   paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                                                   "_RNAMagnet_KEGG_Result_", time_point, ".xlsx"),
                row.names = FALSE, sheetName = paste0("KEGG_Results"))
  }
  
  
  ### DE between two directions in Stroma
  Seurat_Object = Combined_Adult_Seurat_Obj
  target_col = "HSPC"
  comp1 = c("STHSC", "MPP2", "MPP3", "MPP4")
  comp2 = "Stroma"
  time_point = "Adult"
  dim_method = "UMAP"
  result_dir=outputDir2
  
  ### set output directory
  result_dir <- paste0(result_dir, paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                       "_RNAMagnet_Result/", time_point, "/")
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
  
  ### set group info to the metadata
  Seurat_Object@meta.data$Group <- paste0(Seurat_Object@meta.data[,target_col], "_", Seurat_Object@meta.data$Development)
  
  ### all the comps
  comps <- union(paste0(comp1, "_", time_point), paste0(comp2, "_", time_point))
  
  ### set the ident of the object with the specified info
  Seurat_Object <- SetIdent(object = Seurat_Object,
                            cells = rownames(Seurat_Object@meta.data),
                            value = Seurat_Object@meta.data$Group)
  
  ### only keep the specified cells
  Seurat_Object <- subset(Seurat_Object, idents=comps)
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Object)), rownames(Seurat_Object@meta.data)))
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Object@meta.data <- Seurat_Object@meta.data[colnames(Seurat_Object@assays$RNA@counts),]
  
  ### preprocessing
  Seurat_Object <- FindVariableFeatures(Seurat_Object)
  Seurat_Object <- ScaleData(Seurat_Object)
  
  ### run PCA/UMAP
  if(dim_method == "PCA") {
    Seurat_Object <- RunPCA(Seurat_Object, npcs = 15)
    dim_map <- Embeddings(Seurat_Object, reduction = "pca")[rownames(Seurat_Object@meta.data),]
  } else if(dim_method == "UMAP") {
    Seurat_Object <- RunUMAP(Seurat_Object, dims = 1:15)
    dim_map <- Embeddings(Seurat_Object, reduction = "umap")[rownames(Seurat_Object@meta.data),]
  } else {
    stop("ERROR: dim_method not PCA nor UMAP.")
  }
  
  ### run RNAMagnet with anchors
  result <- RNAMagnetAnchors(Seurat_Object,
                             anchors = unique(Seurat_Object@meta.data$Group))
  
  ### add RNAMagnet info to the seurat object
  Seurat_Object@meta.data$direction <- as.character(result[rownames(Seurat_Object@meta.data),"direction"])
  Seurat_Object@meta.data$adhesiveness <- as.numeric(result[rownames(Seurat_Object@meta.data),"adhesiveness"])
  Seurat_Object@meta.data$specificity <- as.numeric(sapply(rownames(Seurat_Object@meta.data), function(x) {
    result[x, result[x,"direction"]]
  }))
  
  ### check the order
  print(identical(rownames(Seurat_Object@meta.data), rownames(dim_map)))
  
  ### set the ident of the object with the Cell Type
  Seurat_Object <- SetIdent(object = Seurat_Object,
                            cells = rownames(Seurat_Object@meta.data),
                            value = Seurat_Object@meta.data$Cell_Type)
  
  ### only keep the specified cells
  Seurat_Object <- subset(Seurat_Object, idents="Stroma")
  
  ### there are many groups in Heme, so unify them
  Seurat_Object@meta.data$direction[which(Seurat_Object@meta.data$direction != "Stroma_Adult")] <- "Other_Adult"
  
  ### set the ident of the object with the Cell Type
  Seurat_Object <- SetIdent(object = Seurat_Object,
                            cells = rownames(Seurat_Object@meta.data),
                            value = Seurat_Object@meta.data$direction)
  
  ### get DE genes between two groups
  de_result <- FindMarkers(object = Seurat_Object, test.use = "wilcox",
                           ident.1 = "Other_Adult", ident.2 = "Stroma_Adult",
                           logfc.threshold = 0, min.pct = 0.1)
  
  ### order the DE reuslt
  de_result <- de_result[order(de_result$p_val_adj),]
  
  ### data frame
  de_result <- data.frame(Gene_Symbol=rownames(de_result),
                          de_result,
                          stringsAsFactors = FALSE, check.names = FALSE)
  
  ### write out the DE result
  write.xlsx2(de_result,
              file = paste0(result_dir,
                            paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                            "_RNAMagnet_DE_Result_", time_point, ".xlsx"),
              sheetName = "DE_Result",
              row.names = FALSE)
  
  ### pathway analysis
  pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Mm.eg.db,
                                                            de_result[which(de_result$p_val_adj < 1e-100),
                                                                      "Gene_Symbol"],
                                                            "ENTREZID", "SYMBOL"),
                                          org = "mouse", database = "GO",
                                          title = paste0(paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                                                         "_RNAMagnet_Pathway_Result"),
                                          displayNum = 30, imgPrint = TRUE,
                                          dir = paste0(result_dir))
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Mm.eg.db,
                                                              de_result[which(de_result$p_val_adj < 1e-100),
                                                                        "Gene_Symbol"],
                                                              "ENTREZID", "SYMBOL"),
                                            org = "mouse", database = "KEGG",
                                            title = paste0(paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                                                           "_RNAMagnet_Pathway_Result"),
                                            displayNum = 30, imgPrint = TRUE,
                                            dir = paste0(result_dir))
  if(!is.null(pathway_result_GO) && nrow(pathway_result_GO) > 0) {
    write.xlsx2(pathway_result_GO, file = paste0(result_dir,
                                                 paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                                                 "_RNAMagnet_GO_Result_", time_point, ".xlsx"),
                row.names = FALSE, sheetName = paste0("GO_Results"))
  }
  if(!is.null(pathway_result_KEGG) && nrow(pathway_result_KEGG) > 0) {
    write.xlsx2(pathway_result_KEGG, file = paste0(result_dir,
                                                   paste(comp1, collapse = "_"), "_vs_", paste(comp2, collapse = "_"),
                                                   "_RNAMagnet_KEGG_Result_", time_point, ".xlsx"),
                row.names = FALSE, sheetName = paste0("KEGG_Results"))
  }
  
}

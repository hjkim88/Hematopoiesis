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
#                                     Robj2_path="./data/Combined_BM_Seurat_Obj.RDATA",
#                                     outputDir="./results/Additional_Oct2020/")
###

additional_analysis <- function(Robj1_path="./data/Combined_Seurat_Obj.RDATA",
                                Robj2_path="./data/Combined_BM_Seurat_Obj.RDATA",
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
  if(!require(slingshot, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("slingshot")
    require(slingshot, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
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
  DimPlot(subset_Seurat_Obj, reduction = "pca", group.by = "Tissue", pt.size = 1.5) +
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
  Seurat_Obj <- subset(Combined_Seurat_Obj, idents=type)
  
  ### order the meta data by developmental time
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[order(factor(Seurat_Obj@meta.data$Development,
                                                            levels = time_points)),]
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@assays$RNA@counts <- Seurat_Obj@assays$RNA@counts[,rownames(Seurat_Obj@meta.data)]
  
  ### run PCA
  Seurat_Obj <- RunPCA(Seurat_Obj, npcs = 10)
  pca_map <- Embeddings(Seurat_Obj, reduction = "pca")[rownames(Seurat_Obj@meta.data),1:10]
  
  ### get slingshot object
  slingshot_obj <- slingshot(pca_map,
                             clusterLabels = Seurat_Obj@meta.data$Development, 
                             reducedDim = "PCA")
  
  
  
  
}

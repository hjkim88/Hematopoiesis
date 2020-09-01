###
#   File name : RNAMagnet_Analysis.R
#   Author    : Hyunjin Kim
#   Date      : Aug 24, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : 1. For Heme objects E16, E18, P0, and Adult, run RNA-Magnet for physical interactions of
#                  each cluster within the object and against the stroma object of the same timepoint.
#                  Also, run RNA-Magnet for physical interactions of each cluster within each stroma object.
#                  * What Trent would like to see is:
#                   a) A tSNE-plot showing the "direction", "adhesiveness", and Specificity Score of the cells
#                      in each Heme object for all anchor populations (Heme and Stroma) at the same timepoint.
#                   b) The same data in a heat map.
#               2. For Heme objects E16, E18, P0, and Adult, run RNA-Magnet to infer the putative signaling
#                  interactions between each cluster within the object and between the stroma objects of the same
#                  timepoint. Also, run RNA-Magnet for putative signaling interactions bewteen cluster within
#                  each stroma object.
#                  * What Trent would like to see here is:
#                   a) A summary of the "specificity" score at the cluster level i) within each Heme object,
#                      ii) between the Heme object and the Stroma object at the same timepoint, and
#                      iii) within each Stroma object, depicted as a signaling network.
#                   b) A list of the molecules within each cluster mediating these interactions.
#                   c) A tSNE-Plot showing the "specificity" of the cells in each Heme object for
#                      i) each cluster within the object and ii) each cluster within the stroma object of the
#                      same timepoint, as well as a tSNE-Plot showing this same value in each Stroma object
#                      for each cluster within the object.
#   
#   Instruction
#               1. Source("RNAMagnet_Analysis.R")
#               2. Run the function "rna_magnet_trent" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_RNAMagnet_Analysis.R/RNAMagnet_Analysis.R")
#               > rna_magnet_trent(Robj_path="C:/Users/hkim8/SJ/HSPC Analyses/Stroma and Heme/",
#                                  Seurat_RObj_path="./data/Combined_BM_Seurat_Obj.RDATA",
#                                  outputDir="./results/RNA-Magnet/")
###

rna_magnet_trent <- function(Robj_path="C:/Users/hkim8/SJ/HSPC Analyses/Stroma and Heme/",
                             Seurat_RObj_path="./data/Combined_BM_Seurat_Obj.RDATA",
                             outputDir="./results/RNA-Magnet/") {
  
  ### load libraries
  options(java.parameters = "-Xmx10240m")
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(remotes, quietly = TRUE)) {
    install.packages("remotes")
    require(remotes, quietly = TRUE)
  }
  if(!require(RNAMagnet, quietly = TRUE)) {
    remotes::install_github("veltenlab/rnamagnet")
    require(RNAMagnet, quietly = TRUE)
  }
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(reticulate, quietly = TRUE)) {
    install.packages("reticulate")
    require(reticulate, quietly = TRUE)
  }
  if(!require(Rmagic, quietly = TRUE)) {
    install.packages("Rmagic")
    require(Rmagic, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(scales, quietly = TRUE)) {
    install.packages("scales")
    library(scales, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  if(!require(gplots, quietly = TRUE)) {
    install.packages("gplots")
    library(gplots, quietly = TRUE)
  }
  if(!require(RColorBrewer, quietly = TRUE)) {
    install.packages("RColorBrewer")
    require(RColorBrewer, quietly = TRUE)
  }
  
  ### see python environment since RNAMagnet uses the python module 'magic'
  # conda_create("r-reticulate")
  # conda_install("r-reticulate", "python-magic")
  use_condaenv("r-reticulate")
  py_config()
  py_module_available("magic")
  
  ### get Robj file list
  f_list <- list.files(Robj_path, pattern = ".Robj$")
  
  ### determine necessary variables
  time_points <- c("E16", "E18", "P0", "ADULT")
  cell_type <- c("Heme", "Stroma")
  possible_names <- as.vector(sapply(cell_type, function(x) paste0(x, time_points)))
  possible_names_mat <- data.frame(Names=possible_names,
                                   Type=as.vector(sapply(cell_type, function(x) rep(x, length(time_points)))),
                                   Time=rep(time_points, length(cell_type)),
                                   stringsAsFactors = FALSE, check.names = FALSE)
  row.names(possible_names_mat) <- possible_names
  
  ### keep the required files only
  f_list <- f_list[which(f_list %in% paste0(possible_names, "_regress.Robj"))]
  
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
  
  ### for each of the Robj file, run the RNAMagnet
  Combined_Seurat_Obj <- NULL
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
    
    ### active assay = "RNA"
    Seurat_Obj@active.assay <- "RNA"
    
    ### set output directory
    outputDir2 <- paste0(outputDir, f, "/")
    dir.create(outputDir2, recursive = TRUE, showWarnings = FALSE)
    
    ### annotate the file name to the R object
    tissue_name <- strsplit(f, split = "_", fixed = TRUE)[[1]][1]
    Seurat_Obj@meta.data$Tissue <- tissue_name
    
    ### annotate the development
    Seurat_Obj@meta.data$Development <- possible_names_mat[tissue_name, "Time"]
    
    ### annotate the Cell type
    Seurat_Obj@meta.data$Cell_Type <- possible_names_mat[tissue_name, "Type"]
    
    ### combine the Seurat R objects
    if(is.null(Combined_Seurat_Obj)) {
      Combined_Seurat_Obj <- Seurat_Obj
    } else {
      Combined_Seurat_Obj <- merge(x = Combined_Seurat_Obj, y = Seurat_Obj,
                                   project = Seurat_Obj@project.name)
    }
    
    ### run only if there are enough number of cells in the R object
    if(nrow(Seurat_Obj@meta.data) > 1 && identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data))) {
      
      ### clustering on the R object
      Seurat_Obj <- FindNeighbors(Seurat_Obj, dims = 1:5, k.param = 5)
      Seurat_Obj <- FindClusters(Seurat_Obj, resolution = 0.4)
      Seurat_Obj@meta.data$new_clusts <- Idents(Seurat_Obj)
      
      #
      ### run RNAMagnet
      #
      ### get PCA matrix
      pca_map <- Embeddings(Seurat_Obj, reduction = "pca")[rownames(Seurat_Obj@meta.data),1:5]
      
      ### set the ident of the seurat object with the cluster info
      Seurat_Obj <- SetIdent(object = Seurat_Obj,
                             cells = rownames(Seurat_Obj@meta.data),
                             value = Seurat_Obj@meta.data$new_clusts)
      
      ### run RNAMagnet
      result <- RNAMagnetAnchors(Seurat_Obj,
                                 anchors = levels(Seurat_Obj@meta.data$new_clusts))
      
      ### make file name
      obj <- strsplit(f, split = "_", fixed = TRUE)[[1]][1]
      
      ### write the result as an Excel file
      write.xlsx2(data.frame(Cell=rownames(result), result,
                             stringsAsFactors = FALSE, check.names = FALSE),
                  file = paste0(outputDir2, "RNAMagnet_Result_", obj, ".xlsx"),
                  sheetName = paste0("RNAMagnet_", obj),
                  row.names = FALSE)
      
      ### add RNAMagnet info to the seurat object
      Seurat_Obj@meta.data$direction <- as.character(result[rownames(Seurat_Obj@meta.data),"direction"])
      Seurat_Obj@meta.data$adhesiveness <- as.numeric(result[rownames(Seurat_Obj@meta.data),"adhesiveness"])
      Seurat_Obj@meta.data$specificity <- as.numeric(sapply(rownames(Seurat_Obj@meta.data), function(x) {
        result[x, paste0("X", result[x,"direction"])]
      }))
      
      ### check the order
      print(identical(rownames(Seurat_Obj@meta.data), rownames(pca_map)))
      
      ### make a data frame for ggplot
      plot_df <- data.frame(X=pca_map[rownames(Seurat_Obj@meta.data),"PC_1"],
                            Y=pca_map[rownames(Seurat_Obj@meta.data),"PC_2"],
                            group_color = Seurat_Obj@meta.data$direction,
                            group_alpha = Seurat_Obj@meta.data$adhesiveness,
                            cluster_color = Seurat_Obj@meta.data$new_clusts,
                            specificity = Seurat_Obj@meta.data$specificity,
                            stringsAsFactors = FALSE, check.names = FALSE)
      
      ### factorize the direction columns
      plot_df$group_color <- factor(plot_df$group_color, levels = levels(plot_df$cluster_color)) 
      
      ### get colors for the clustering result
      cell_colors_clust <- cell_pal(levels(Seurat_Obj@meta.data$new_clusts), hue_pal())
      
      ### scatter plot
      p <- list()
      
      ### original PCA
      p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
        geom_point(aes_string(col="cluster_color"), size=2) +
        xlab("PC1") + ylab("PC2") +
        labs(col="Cluster") +
        ggtitle("PCA with Cluster Info") +
        theme_classic(base_size = 16) +
        scale_color_manual(values = cell_colors_clust, labels = names(cell_colors_clust))
      p[[1]] <- LabelClusters(plot = p[[1]], id = "cluster_color")
      
      ### draw a scatter plot with the adhesiveness info
      p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
        geom_point(aes_string(col="group_color", alpha="group_alpha"), size=2) +
        xlab("PC1") + ylab("PC2") +
        labs(col="Direction", alpha="Adhesiveness") +
        ggtitle("PCA with Direction & Adhesiveness") +
        theme_classic(base_size = 16) +
        scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$group_color)[order(unique(plot_df$group_color))])],
                           labels = names(cell_colors_clust[as.character(unique(plot_df$group_color)[order(unique(plot_df$group_color))])]))
      p[[2]] <- LabelClusters(plot = p[[2]], id = "cluster_color")
      
      ### save the plots
      g <- arrangeGrob(grobs = p,
                       nrow = 2,
                       ncol = 1,
                       top = "")
      ggsave(file = paste0(outputDir2, "PCA_RNAMagnet_Result_AD_", obj, ".png"), g, width = 20, height = 12, dpi = 300)
      
      ### draw a scatter plot with the specificity info
      p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
        geom_point(aes_string(col="group_color", alpha="specificity"), size=2) +
        xlab("PC1") + ylab("PC2") +
        labs(col="Direction", alpha="Specificity Score") +
        ggtitle("PCA with Direction & Specificity Score") +
        theme_classic(base_size = 16) +
        scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$group_color)[order(unique(plot_df$group_color))])],
                           labels = names(cell_colors_clust[as.character(unique(plot_df$group_color)[order(unique(plot_df$group_color))])]))
      p[[2]] <- LabelClusters(plot = p[[2]], id = "cluster_color")
      
      ### save the plots
      g <- arrangeGrob(grobs = p,
                       nrow = 2,
                       ncol = 1,
                       top = "")
      ggsave(file = paste0(outputDir2, "PCA_RNAMagnet_Result_SP_", obj, ".png"), g, width = 20, height = 12, dpi = 300)
      
      #
      ### heatmap with specificity scores - col side color bars with direction & adhesiveness
      #
      
      ### get a matrix for the heatmap
      heatmap_mat <- data.frame(result[,paste0("X", levels(Seurat_Obj@meta.data$new_clusts))],
                                stringsAsFactors = FALSE, check.names = FALSE)
      
      ### scale the data
      heatmap_mat_scaled <- scale_h(heatmap_mat, type = "row")
      
      ### set rowside colors1 (direction)
      ### and it can also be used as rowside colors
      uniqueV <- levels(Seurat_Obj@meta.data$new_clusts)
      # colors1 <- colorRampPalette(brewer.pal(9,"Blues"))(length(uniqueV))
      colors1 <- rainbow(length(uniqueV))
      names(colors1) <- uniqueV
      
      ### set rowside colors2 (adhesiveness)
      uniqueV <- unique(result$adhesiveness)[order(unique(result$adhesiveness))]
      # colors2 <- colorRampPalette(brewer.pal(9,"Blues"))(length(uniqueV))
      colors2 <- gray.colors(length(uniqueV))
      names(colors2) <- uniqueV
      
      ### hierarchical clustering functions
      dist.spear <- function(x) as.dist(1-cor(t(x), method = "spearman"))
      hclust.ave <- function(x) hclust(x, method="average")
      
      ### heatmap
      rowSide <- t(cbind(colors2[as.character(result$adhesiveness)],
                         colors1[as.character(Seurat_Obj@meta.data[rownames(result),"new_clusts"])],
                         colors1[as.character(result$direction)]))
      rownames(rowSide) <- c("Adhesiveness", "Cluster", "Direction")
      png(paste0(outputDir2, "Heatmap_RNAMagnet_Result_", obj, ".png"), width = 2400, height = 1400, res = 120)
      par(oma=c(4,0,0,0))
      heatmap.3(as.matrix(heatmap_mat_scaled), main = paste0("RNAMagnet_Specificity_", obj),
                xlab = "Clusters", ylab = "", col=greenred(100),
                scale="none", key=T, keysize=0.8, density.info="density",
                dendrogram = "none", trace = "none",
                labRow = FALSE, labCol = levels(Seurat_Obj@meta.data$new_clusts),
                Rowv = TRUE, Colv = FALSE,
                distfun=dist.spear, hclustfun=hclust.ave,
                ColSideColors = cbind(colors1[levels(Seurat_Obj@meta.data$new_clusts)]),
                RowSideColors = rowSide,
                cexRow = 3, cexCol = 3, na.rm = TRUE)
      ### legend
      lgd = rep(NA, 19)
      lgd[c(1,10,19)] = c(signif(max(as.numeric(result$adhesiveness)), 3),
                        signif(mean(as.numeric(result$adhesiveness)), 3),
                        signif(min(as.numeric(result$adhesiveness)), 3))
      legend("left",
             title = "Adhesiveness\n\n",
             legend = lgd,
             fill = gray.colors(19),
             border = NA,
             bty = "o",
             x.intersp = 1,
             y.intersp = 0.5,
             cex = 1)
      dev.off()
      
    }
    
  }
  
  ### set the ident of the object with the Development
  Combined_Seurat_Obj <- SetIdent(object = Combined_Seurat_Obj,
                                  cells = rownames(Combined_Seurat_Obj@meta.data),
                                  value = Combined_Seurat_Obj@meta.data$Development)
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Combined_Seurat_Obj)), rownames(Combined_Seurat_Obj@meta.data)))
  
  ### for each time run RNAMagnet between Heme and Stroma
  for(tp in time_points) {
    
    ### new output directory
    outputDir2 <- paste0(outputDir, tp, "/")
    dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
    
    ### split the Seurat obj based on Cell Type info
    Seurat_Obj <- subset(Combined_Seurat_Obj, idents=tp)
    
    ### order the meta data by cell type
    Seurat_Obj@meta.data <- Seurat_Obj@meta.data[order(Seurat_Obj@meta.data$Cell_Type),]
    
    ### rownames in the meta.data should be in the same order as colnames in the counts
    Seurat_Obj@assays$RNA@counts <- Seurat_Obj@assays$RNA@counts[,rownames(Seurat_Obj@meta.data)]
    
    ### run PCA
    Seurat_Obj <- RunPCA(Seurat_Obj, npcs = 15)
    pca_map <- Embeddings(Seurat_Obj, reduction = "pca")[rownames(Seurat_Obj@meta.data),1:15]
    
    ### clustering on the R object
    Seurat_Obj <- FindNeighbors(Seurat_Obj, dims = 1:5, k.param = 5)
    Seurat_Obj <- FindClusters(Seurat_Obj, resolution = 0.4)
    Seurat_Obj@meta.data$new_clusts <- Idents(Seurat_Obj)
    
    ### set the ident of the seurat object with the cluster info
    Seurat_Obj <- SetIdent(object = Seurat_Obj,
                           cells = rownames(Seurat_Obj@meta.data),
                           value = Seurat_Obj@meta.data$new_clusts)
    
    ### run RNAMagnet
    result <- RNAMagnetAnchors(Seurat_Obj,
                               anchors = levels(Seurat_Obj@meta.data$new_clusts))
    
    ### write the result as an Excel file
    write.xlsx2(data.frame(Cell=rownames(result), result,
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(outputDir2, "RNAMagnet_Result_", tp, ".xlsx"),
                sheetName = paste0("RNAMagnet_", tp),
                row.names = FALSE)
    
    ### add RNAMagnet info to the seurat object
    Seurat_Obj@meta.data$direction <- as.character(result[rownames(Seurat_Obj@meta.data),"direction"])
    Seurat_Obj@meta.data$adhesiveness <- as.numeric(result[rownames(Seurat_Obj@meta.data),"adhesiveness"])
    Seurat_Obj@meta.data$specificity <- as.numeric(sapply(rownames(Seurat_Obj@meta.data), function(x) {
      result[x, paste0("X", result[x,"direction"])]
    }))
    
    ### check the order
    print(identical(rownames(Seurat_Obj@meta.data), rownames(pca_map)))
    
    ### make a data frame for ggplot
    plot_df <- data.frame(X=pca_map[rownames(Seurat_Obj@meta.data),"PC_1"],
                          Y=pca_map[rownames(Seurat_Obj@meta.data),"PC_2"],
                          group_color = Seurat_Obj@meta.data$direction,
                          group_alpha = Seurat_Obj@meta.data$adhesiveness,
                          cluster_color = Seurat_Obj@meta.data$new_clusts,
                          specificity = Seurat_Obj@meta.data$specificity,
                          cell_type = Seurat_Obj@meta.data$Cell_Type,
                          stringsAsFactors = FALSE, check.names = FALSE)
    
    ### factorize the direction columns
    plot_df$group_color <- factor(plot_df$group_color, levels = levels(plot_df$cluster_color)) 
    
    ### get colors for the clustering result
    cell_colors_clust <- cell_pal(levels(Seurat_Obj@meta.data$new_clusts), hue_pal())
    
    ### scatter plot
    p <- list()
    
    ### original PCA
    p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="cell_type"), size=2) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Cluster") +
      ggtitle("PCA with Cell Type") +
      theme_classic(base_size = 16)
    p[[1]] <- LabelClusters(plot = p[[1]], id = "cluster_color")
    
    ### draw a scatter plot with the adhesiveness info
    p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="group_color", alpha="group_alpha"), size=2) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Direction", alpha="Adhesiveness") +
      ggtitle("PCA with Direction & Adhesiveness") +
      theme_classic(base_size = 16) +
      scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$group_color)[order(unique(plot_df$group_color))])],
                         labels = names(cell_colors_clust[as.character(unique(plot_df$group_color)[order(unique(plot_df$group_color))])]))
    p[[2]] <- LabelClusters(plot = p[[2]], id = "cluster_color")
    
    ### save the plots
    g <- arrangeGrob(grobs = p,
                     nrow = 2,
                     ncol = 1,
                     top = "")
    ggsave(file = paste0(outputDir2, "PCA_RNAMagnet_Result_AD_", tp, ".png"), g, width = 20, height = 12, dpi = 300)
    
    ### draw a scatter plot with the specificity info
    p[[2]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="group_color", alpha="specificity"), size=2) +
      xlab("PC1") + ylab("PC2") +
      labs(col="Direction", alpha="Specificity Score") +
      ggtitle("PCA with Direction & Specificity Score") +
      theme_classic(base_size = 16) +
      scale_color_manual(values = cell_colors_clust[as.character(unique(plot_df$group_color)[order(unique(plot_df$group_color))])],
                         labels = names(cell_colors_clust[as.character(unique(plot_df$group_color)[order(unique(plot_df$group_color))])]))
    p[[2]] <- LabelClusters(plot = p[[2]], id = "cluster_color")
    
    ### save the plots
    g <- arrangeGrob(grobs = p,
                     nrow = 2,
                     ncol = 1,
                     top = "")
    ggsave(file = paste0(outputDir2, "PCA_RNAMagnet_Result_SP_", tp, ".png"), g, width = 20, height = 12, dpi = 300)
    
    #
    ### heatmap with specificity scores - col side color bars with direction & adhesiveness
    #
    
    ### get a matrix for the heatmap
    heatmap_mat <- data.frame(result[,paste0("X", levels(Seurat_Obj@meta.data$new_clusts))],
                              stringsAsFactors = FALSE, check.names = FALSE)
    
    ### scale the data
    heatmap_mat_scaled <- scale_h(heatmap_mat, type = "row")
    
    ### set rowside colors1 (direction)
    ### and it can also be used as rowside colors
    uniqueV <- levels(Seurat_Obj@meta.data$new_clusts)
    # colors1 <- colorRampPalette(brewer.pal(9,"Blues"))(length(uniqueV))
    colors1 <- rainbow(length(uniqueV))
    names(colors1) <- uniqueV
    
    ### set rowside colors2 (adhesiveness)
    uniqueV <- unique(result$adhesiveness)[order(unique(result$adhesiveness))]
    # colors2 <- colorRampPalette(brewer.pal(9,"Blues"))(length(uniqueV))
    colors2 <- gray.colors(length(uniqueV))
    names(colors2) <- uniqueV
    
    ### set rowside colors3 (cell type)
    uniqueV <- unique(Seurat_Obj@meta.data$Cell_Type)
    colors3 <- hue_pal()(2)
    names(colors3) <- uniqueV
    
    ### hierarchical clustering functions
    dist.spear <- function(x) as.dist(1-cor(t(x), method = "spearman"))
    hclust.ave <- function(x) hclust(x, method="average")
    
    ### heatmap
    rowSide <- t(cbind(colors2[as.character(result$adhesiveness)],
                       colors1[as.character(Seurat_Obj@meta.data[rownames(result),"new_clusts"])],
                       colors1[as.character(result$direction)],
                       colors3[as.character(Seurat_Obj@meta.data[rownames(result),"Cell_Type"])]))
    rownames(rowSide) <- c("Adhesiveness", "Cluster", "Direction", "Cell Type")
    png(paste0(outputDir2, "Heatmap_RNAMagnet_Result_", tp, ".png"), width = 3000, height = 1200, res = 120)
    par(oma=c(4,0,0,0))
    heatmap.3(as.matrix(heatmap_mat_scaled), main = paste0("RNAMagnet_Specificity_", tp),
              xlab = "Clusters", ylab = "", col=greenred(100),
              scale="none", key=T, keysize=0.8, density.info="density",
              dendrogram = "none", trace = "none",
              labRow = FALSE, labCol = levels(Seurat_Obj@meta.data$new_clusts),
              Rowv = TRUE, Colv = FALSE,
              distfun=dist.spear, hclustfun=hclust.ave,
              ColSideColors = cbind(colors1[levels(Seurat_Obj@meta.data$new_clusts)]),
              RowSideColors = rowSide,
              cexRow = 3, cexCol = 3, na.rm = TRUE)
    ### legend
    lgd = rep(NA, 19)
    lgd[c(1,10,19)] = c(signif(max(as.numeric(result$adhesiveness)), 3),
                        signif(mean(as.numeric(result$adhesiveness)), 3),
                        signif(min(as.numeric(result$adhesiveness)), 3))
    legend("left",
           title = "Adhesiveness\n\n",
           legend = lgd,
           fill = gray.colors(19),
           border = NA,
           bty = "o",
           x.intersp = 1,
           y.intersp = 0.5,
           cex = 1)
    legend("bottomleft",
           title = "Cell Type",
           legend = names(colors3),
           fill = colors3,
           border = NA,
           bty = "o",
           x.intersp = 1,
           y.intersp = 1,
           cex = 1)
    dev.off()
    
  }
  
}

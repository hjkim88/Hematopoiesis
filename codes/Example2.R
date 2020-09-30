### THIS IS AN EXAMPLE CODE FOR PC1 GENE - HEATMAP & PATHWAY ANALYSIS

### load libraries
if(!require(Seurat, quietly = TRUE)) {
  install.packages("Seurat")
  require(Seurat, quietly = TRUE)
}
if(!require(ggplot2, quietly = TRUE)) {
  install.packages("ggplot2")
  require(ggplot2, quietly = TRUE)
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
#
### MPP2
#
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

### PC1 density plot
plot(density(pca_contb[,"PC_1"]))

### get genes that contributed to the PC1 the most
contb_threshold <- 0.1
important_genes <- rownames(pca_contb)[which(pca_contb[,"PC_1"] > contb_threshold)]

### get genes that contributed to the PC1 the least
garbage_threshold <- 1e-04
garbage_genes <- rownames(pca_contb)[which(pca_contb[,"PC_1"] < garbage_threshold)]


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

#
### a heatmap with the important genes - 322 genes x 428 cells
#

### get a matrix for the heatmap
heatmap_mat <- data.frame(Seurat_Obj@assays$RNA@counts[important_genes,], check.names = FALSE)

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


#
### a heatmap with the garbage genes - 167 genes x 428 cells
#

### get a matrix for the heatmap
heatmap_mat <- data.frame(Seurat_Obj@assays$RNA@counts[garbage_genes,], check.names = FALSE)

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

### heatmap
png(paste0(outputDir2, "PC1_Genes_Heatmap_GARBAGE.png"), width = 2000, height = 1000)
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

#
### msigdb pacakge + fgsea
### signature: PC1 contribution
#

### signature preparation
signat <- pca_contb[important_genes,"PC_1"]
names(signat) <- important_genes

### db preparation
# MSIGDB
m_df <- msigdbr(species = "Mus musculus")
m_list <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

### run GSEA
GSEA_result <- run_gsea(gene_list = m_list, signature = list(signat), printPlot = FALSE)
GSEA_result <- GSEA_result[order(GSEA_result$pval),]

### only get pathways that have pval < 0.05
pathways <- GSEA_result$pathway[intersect(which(GSEA_result$pval < 0.05),
                                          which(GSEA_result$size > 10))]

### run GSEA again with the significant result - plot printing
GSEA_result2 <- run_gsea(gene_list = m_list[pathways], signature = list(signat),
                         printPlot = TRUE, printPath = "./results/GSEA/")
GSEA_result2 <- GSEA_result2[order(GSEA_result2$padj, GSEA_result2$pval),]

### write out the result file
write.xlsx2(GSEA_result2, file = paste0("./results/GSEA/MPP2_GSEA_msigdb_result.xlsx"),
            sheetName = "GSEA_Result", row.names = FALSE)






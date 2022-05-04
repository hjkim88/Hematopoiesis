### A sript for checking # cell and read depth of the new hematopoietic samples

### load libraries
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

### sample location
sample_dir <- "Z:/ResearchHome/SharedResources/Immunoinformatics/common/JCC/JCC395/"

### get folders in there
sample_folders <- list.dirs(path = sample_dir, full.names = FALSE, recursive = FALSE)

### load metric files into R memory space
metric_info <- vector("list", length = length(sample_folders))
names(metric_info) <- sample_folders
metric_info2 <- NULL
for(sp in names(metric_info)) {
  writeLines(paste(sp))
  metric_info[[sp]] <- t(read.csv(file = paste0(sample_dir, sp, "/metrics_summary.csv"),
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
                      Antibody_Number_of_Reads=as.numeric(gsub(",", "", metric_info2["Antibody: Number of Reads",])),
                      Antibody_Reads_in_Cells=as.numeric(gsub("%", "", metric_info2["Antibody: Antibody Reads in Cells",])),
                      Antibody_Mean_Reads_per_Cell=as.numeric(gsub(",", "", metric_info2["Antibody: Mean Reads per Cell",])),
                      Antibody_Median_UMIs_per_Cell=as.numeric(gsub(",", "", metric_info2["Antibody: Median UMIs per Cell (summed over all recognized antibody barcodes)",])),
                      stringsAsFactors = FALSE, check.names = FALSE)

p <- vector("list", length = ncol(plot_df)-1)
names(p) <- colnames(plot_df)[-1]

for(col in names(p)) {
  p[[col]] <- ggplot(data=plot_df, aes_string(x="Sample", y=col, label=col)) +
    geom_bar(stat = "identity", fill = "gray60") +
    ggtitle("") +
    geom_text(hjust=0.6, color="black", size=3) +
    coord_flip() +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title = element_text(hjust = 0.5, vjust = 0.5, size = 24),
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16),
          axis.ticks = element_blank())
}


g <- arrangeGrob(grobs = p,
                 nrow = 4,
                 ncol = 3)
ggsave(file = paste0("./results/New_Samples/QC_Results.pdf"), g, width = 25, height = 18)



###
#   File name : Pseudotime_Analysis.R
#   Author    : Hyunjin Kim
#   Date      : Aug 8, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Pseudotime Analyses on the clusters within each R object.
#               Then Combine all the R object and do the same analysis across the development
#               (E16LTHSC, E18LTHSC, P0LTHSC, and ADULTLTHSC, etc.)
#   
#   Instruction
#               1. Source("Pseudotime_Analysis.R")
#               2. Run the function "go_analysis_trent" - specify the input directory and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Pseudotime_Analysis.R/Pseudotime_Analysis.R")
#               > go_analysis_trent(Robj_path="C:/Users/hkim8/SJ/HSPC Analyses/HSPC Subsets/",
#                                   outputDir="./results/")
###

go_analysis_trent <- function(Robj_path="C:/Users/hkim8/SJ/HSPC Analyses/HSPC Subsets/",
                              outputDir="./results/Pseudotime/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  
  ### get Robj file list
  f_list <- list.files(Robj_path, pattern = ".Robj$")
  
  ### for each of the Robj file, run the pseudotime analysis and combine the R objects for later use
  for(f in f_list) {
    
    
    
    
  }
  
  
  
  
}

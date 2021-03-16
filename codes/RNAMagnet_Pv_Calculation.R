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
  
  
  
  
  AdultStromavLTHSC_NA_Anchors <- RNAMagnetAnchors(AdultStromavLTHSC_NA, c("LTHSCs",'Stroma'), return = "summary", neighborhood.distance = 0.7, neighborhood.gradient = 3, .k = 10, .x0 = 0.5, .minExpression = 0, .minMolecules = 1, .version = "1.0.0", .cellularCompartment = c("Membrane","ECM","Both"), .manualAnnotation = "Correct" )
  
  
  
}





# Testscript for using the R package cellPixels for development  ++++++++++
# Author: Kai Budde
# Last changed: 2020/11/24

#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


# Delete everything in the environment
rm(list = ls())
# close all open plots in RStudio
graphics.off()


# Please adapt the following parameters ####################################
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Directory of the images
input_dir <- "test2/"

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Load packages #############################################################

if(!("EBImage" %in% installed.packages())){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("EBImage")
  require("EBImage")
}


#list.of.packages <- c("tiff", "dplyr", "devtools", "xlsx")
list.of.packages <- c("tiff", "devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require(tiff)
# require(dplyr)
require(devtools)
# require(EBImage)
# require(xlsx)
#
# # Install the R package for producing stacks of the images
# if(!("stackImages" %in% installed.packages()[,"Package"])){
#   if(!installed.packages()[,"Version"][installed.packages()[,"Package"] == "stackImages"] == "0.1.4"){
#     devtools::install_github("SFB-ELAINE/stackImages", ref = "v0.1.4")
#   }
# }
# require(stackImages)

# Check package
check()

# Document package
document()

# Load package to use it
load_all()

## FIRST EXAMPLE DIRECTORY -------------------------------------------------




# Data frame with certain parameter values
#number_of_dirs <- length(input_dirs)
#df_results <- data.frame("Directory" = input_dirs,
#                         "threshold_find" = rep(-99, number_of_dirs),
#                         "threshold_connect" = rep(-99, number_of_dirs)
#                         )

# # Obtain all positions in every z-layer and legths of al cilia
# for(i in 1:number_of_dirs){
#   input_dir <- input_dirs[i]
#   output_list <- detectCilia(input_dir = input_dir,
#                              cilium_color = cilium_color,
#                              threshold_by_density_of_cilium_pixels = threshold_by_density_of_cilium_pixels,
#                              vicinity = vicinity,
#                              min_size = min_size,
#                              max_size = max_size,
#                              number_size_factor = number_size_factor,
#                              pixel_size = pixel_size,
#                              slice_distance = slice_distance)
#   if(!is.null(output_list)){
#     df_results$threshold_find[df_results$Directory == input_dir] <-
#       output_list$df_parameterlist$parameterValues[
#         output_list$df_parameterlist$parameterNames == "threshold_find"]
#     df_results$threshold_connect[df_results$Directory == input_dir] <-
#       output_list$df_parameterlist$parameterValues[
#         output_list$df_parameterlist$parameterNames == "threshold_connect"]
#   }
#
# }


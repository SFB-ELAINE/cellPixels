# Testscript for using the R package cellPixels for development  +++++++++++
# Author: Kai Budde
# Last changed: 2020/12/04

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

list.of.packages <- c("tiff", "devtools", "reticulate")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require(tiff)
require(devtools)
require(reticulate)

# Check package
check()

# Document package
document()

# Load package to use it
load_all()

## FIRST EXAMPLE DIRECTORY -------------------------------------------------

df_results <- cellPixels(input_dir = input_dir,
                         nucleus_color = "blue",
                         number_size_factor = 0.2)

# Testscript for using the R package cellPixels for development  +++++++++++
# Author: Kai Budde
# Created:
# Last changed: 2021/10/14

#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


# TODO: Im Namen sollte ersichtlich sein, welche Farbe gesucht und markiert wurde.

# Delete everything in the environment
rm(list = ls())
# close all open plots in RStudio
graphics.off()

# Please adapt the following parameters ####################################
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Directory of the images
#input_dir <- "testleonie/"
#input_dir <- "test4/"
#input_dir <- "E:/bcat/"
#input_dir <- "E:/ROX/"
#input_dir <- "E:/Cav/"
input_dir <- "E:/LSM/"
input_dir <- "E:/test/"
#input_dir <- "testapotome/"

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Load packages #############################################################

if(!("EBImage" %in% installed.packages())){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("EBImage")
}
require("EBImage")

list.of.packages <- c("tiff", "devtools", "reticulate")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  install.packages(new.packages)
}
require(tiff)
require(devtools)
require(reticulate)


# Install readCzi if not already installed (checked with v0.1.9)
devtools::install_github(repo = "https://github.com/SFB-ELAINE/readCzi")

require(readCzi)

# Check package
#check()

# Document package
document()

# Load package to use it
load_all()

## FIRST EXAMPLE DIRECTORY -------------------------------------------------

#BCAT
df_results <- cellPixels(input_dir = input_dir,
                        nucleus_color = "blue",
                        protein_in_nuc_color = "red",
                        protein_in_cytosol_color = "green",
                        number_size_factor = 0.2,
                        add_scale_bar = TRUE)

#ROX
# df_results <- cellPixels(input_dir = input_dir,
#                          nucleus_color = "blue",
#                          protein_in_nuc_color = "green",
#                          protein_in_cytosol_color = "red",
#                          number_size_factor = 0.2)

#Caveolin (Apotome)
# df_results <- cellPixels(input_dir = input_dir,
#                          apotome = TRUE,
#                          apotome_section = FALSE,
#                          nucleus_color = "blue",
#                          protein_in_nuc_color = "none",
#                          protein_in_cytosol_color = "green",
#                          protein_in_membrane_color = "red",
#                          number_size_factor = 0.2)

# df_results <- cellPixels(input_dir = input_dir,
#                          apotome = TRUE,
#                          apotome_section = TRUE,
#                          nucleus_color = "blue",
#                          protein_in_nuc_color = "none",
#                          protein_in_cytosol_color = "green",
#                          protein_in_membrane_color = "red",
#                          number_size_factor = 0.2)


# df_results <- cellPixels(input_dir = input_dir,
#                          nucleus_color = "green",
#                          protein_in_nuc_color = "none",
#                          protein_in_cytosol_color = "none",
#                          number_size_factor = 0.2,
#                          thresh_w_h_nuc = 5,
#                          thresh_offset = 0.0001,
#                          blur_sigma = 0,
#                          use_histogram_equalized = TRUE)

save(df_results, file=paste(input_dir, "output/df_results.Rda", sep=""))

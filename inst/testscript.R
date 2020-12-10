# Testscript for using the R package cellPixels ++++++++++++++++++++++++++++
# Author: Kai Budde
# Last changed: 2020/12/10

# Delete everything in the environment
rm(list = ls())
# close all open plots in RStudio
graphics.off()


# Please adapt the following parameters ####################################
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Directory of the images
input_dir <- "test/"

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Load packages ############################################################

list.of.packages <- c("devtools")
new.packages <- list.of.packages[
  !(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require(devtools)

# Install the R package for producing stacks of the images
devtools::install_github("SFB-ELAINE/cellPixels")
require(cellPixels)

## FIRST EXAMPLE DIRECTORY -------------------------------------------------

df_results <- cellPixels(input_dir = input_dir,
                        nucleus_color = "blue",
                        protein_in_nuc_color = "red",
                        number_size_factor = 0.2)

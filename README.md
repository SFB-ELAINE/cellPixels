# cellPixels

R package for statistics of fluorescence microscopy images.
The package need at least R 4.0.0.
It uses the R packages "reticulate" to load and use python packages (this will install miniconda when first used),
"EBImage" for finding nuclei and contrast enhancement as well as "tiff" for saving the results in the tif format.
We are using the python package "czifile" to directly work with the microscopy image format czi.
Important: Because this pacakge makes use of a python package which is being installed on-the-fly, you need internet connection.


## Code for using the R package

The most recent version of the following code can always be found in
[/inst/testscript.R](https://github.com/SFB-ELAINE/cellPixels/blob/main/inst/testscript.R).

```R
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
```

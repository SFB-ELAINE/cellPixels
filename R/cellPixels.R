#' @title cellPixels
#' @description Main function for counting pixels in different regions
#' @details Input should be tif-format.
#' We are counting the intensities of pixels of different color within and
#' outside of the nuclei. We also add information about the entire image.
#' @aliases cellpixels CellPixels
#' @author Kai Budde
#' @export cellPixels
#' @param input_dir A character (directory that contains all images)

cellPixels <- function(input_dir = NULL,
                       nucleus_color = "blue",
                       number_size_factor = 0.1) {

  # Basics and sourcing functions ------------------------------------------
  .old.options <- options()
  on.exit(options(.old.options))

  options(stringsAsFactors = FALSE, warn=-1)

  # ---------------------------------------------------------------------- #
  # ---------------------- Data acquisition ------------------------------ #
  # ---------------------------------------------------------------------- #

  # Check for NULLs in parameter values ------------------------------------

  # Input directory must be submitted. If not: close function call.
  if(is.null(input_dir)){
    print(paste("Please call the function with an input directory ",
                "which contains fluorescence microscopy images of cells",
                sep=""))
    return()
  }

  # Data input and output --------------------------------------------------

  # Save the file names (tifs) ---------------------------------------------
  file_names <- list.files(path = input_dir)
  file_names <- file_names[grepl("tif", file_names)]
  number_of_tifs <- length(file_names)

  # If there is now tif-file, close function call
  if(number_of_tifs == 0){
    return()
  }

  # Make a new subdirectory inside the input directory
  if(grepl("\\\\", input_dir)){
    input_dir <- gsub("\\$", "", input_dir)
    output_dir <- paste(input_dir, "\\output\\", sep="")
  }else{
    input_dir <- gsub("/$", "", input_dir)
    output_dir <- paste(input_dir, "/output/", sep="")
  }

  dir.create(output_dir, showWarnings = FALSE)

  # Create empty data fram -------------------------------------------------
  df_results <- data.frame(
    "fileName" = file_names,
    "dimension_x" = rep(NA, number_of_tifs),
    "dimension_y" = rep(NA, number_of_tifs),
    "intensity_sum_red_full" = rep(NA, number_of_tifs),
    "intensity_sum_green_full" = rep(NA, number_of_tifs),
    "intensity_sum_blue_full" = rep(NA, number_of_tifs),
    "intensity_sum_red_nucleus_region" = rep(NA, number_of_tifs),
    "intensity_sum_green_nucleus_region" = rep(NA, number_of_tifs),
    "intensity_sum_blue_nucleus_region" = rep(NA, number_of_tifs),
    "intensity_sum_red_without_nucleus_region" = rep(NA, number_of_tifs),
    "intensity_sum_green_without_nucleus_region" = rep(NA, number_of_tifs),
    "intensity_sum_blue_without_nucleus_region" = rep(NA, number_of_tifs))

  # Go through every image in the directory --------------------------------
  for(i in 1:number_of_tifs){

    print(paste("Dealing with file >>", file_names[i], "<<. (It is now ",
                Sys.time(), ".)", sep=""))

    # Get the image name without the ".tif" ending
    image_name_wo_tif <- gsub("\\.tif", "", file_names[i])

    # Get the image path
    image_path <- paste(input_dir, file_names[i], sep="/")

    # Load image
    image_loaded <- tiff::readTIFF(source = image_path, info = FALSE)

    # -------------------------------------------------------------------- #
    # ---------------------- Data manipulation --------------------------- #
    # -------------------------------------------------------------------- #

    # Normalize intensity --------------------------------------------------
    image_normalized <- normalizeIntensity(image = image_loaded)

    # Use Contrast Limited Adaptive Histogram Equalization
    image_histogram_equalization <- EBImage::clahe(x = image_loaded, nx = 4)

    # Find the nuclei ------------------------------------------------------

    # Save only color layer of nuclei
    image_nuclei <- getLayer(image = image_loaded, layer = nucleus_color)
    image_nuclei <- EBImage::Image(image_nuclei)
    #display(image_nuclei)

    # Blur the image
    image_nuclei <- EBImage::gblur(image_nuclei, sigma = 5)

    # Mask the nuclei
    nmask <- EBImage::thresh(image_nuclei, w=15, h=15, offset=0.05)

        # Morphological opening to remove objects smaller than the structuring element
    nmask <- EBImage::opening(nmask, makeBrush(5, shape='disc'))
    # Fill holes
    nmask <- EBImage::fillHull(nmask)

    # Save black-and-white-image as nucleus mask
    nucleus_mask <- nmask

    # Label each connected set of pixels with a distinct ID
    nmask <- EBImage::bwlabel(nmask)

    #display(nmask)

    # Record all nuclei that are at the edges of the image
    left  <- table(nmask[1, 1:dim(nmask)[2]])
    top   <- table(nmask[1:dim(nmask)[1],1])
    right <- table(nmask[dim(nmask)[1], 1:dim(nmask)[2]])
    bottom <- table(nmask[1:dim(nmask)[1],dim(nmask)[2]])

    left <- as.integer(names(left))
    top <- as.integer(names(top))
    right <- as.integer(names(right))
    bottom <- as.integer(names(bottom))

    nuclei_at_borders <- unique(c(left, top, right, bottom))

    # delete the 0 in the list (these are pixels that do not belong to a nucleus)
    nuclei_at_borders <- nuclei_at_borders[nuclei_at_borders != 0]

    # Delete all nuclei at border
    if(length(nuclei_at_borders) > 0){
      for(i in 1:length(nuclei_at_borders)){
        imageData(nmask)[imageData(nmask) == nuclei_at_borders[i]] <- 0
      }
    }

    #display(nmask)

    # Delete all remaining nuclei that are smaller than 5% of the median size
    # object sizes
    # barplot(table(nmask)[-1])

    table_nmask <- table(nmask)
    nuc_min_size <- 0.05*median(table_nmask[-1])

    # remove objects that are smaller than min_nuc_size
    to_be_removed <- as.integer(names(which(table_nmask < nuc_min_size)))
    if(length(to_be_removed) > 0){
      for(i in 1:length(to_be_removed)){
        imageData(nmask)[imageData(nmask) == to_be_removed[i]] <- 0
      }
    }

    # Recount nuclei
    nmask <- EBImage::bwlabel(nmask)
    #display(nmask)

    # Watershed in order to distinct nuclei that are too close to each other
    nmask_watershed <-  EBImage::watershed( distmap(nmask), tolerance = 0.5, ext = 1)
    #display(colorLabels(nmask_watershed), all=TRUE)

    # Count number of cells
    nucNo <- max(nmask_watershed)

    # Include numbers of nuclei
    image_nuclei <- image_loaded
    image_nuclei_numbers <- image_loaded

    table_nmask_watershed <- table(nmask_watershed)

    if(length(table_nmask_watershed[-1]) > 0){
      # remove 0
      nuc_numbers <- as.integer(names(table_nmask_watershed[-1]))

      for(i in 1:length(nuc_numbers)){

        # Find approximate midpoint of every nucleus
        dummy_coordinates <- which(
          imageData(nmask_watershed) == nuc_numbers[i], arr.ind = TRUE)


        pos_x <- round(mean(dummy_coordinates[,1]))
        pos_y <- round(mean(dummy_coordinates[,2]))

        image_nuclei_numbers <- addNumberToImage(image = image_nuclei_numbers,
                                                number = i,
                                                pos_x = pos_x,
                                                pos_y = pos_y,
                                                number_size_factor = number_size_factor,
                                                number_color = "green")
        image_nuclei_numbers <- addNumberToImage(image = image_nuclei_numbers,
                                                number = i,
                                                pos_x = pos_x,
                                                pos_y = pos_y,
                                                number_size_factor = number_size_factor,
                                                number_color = "blue")
      }
    }

    # Add border of nuclei and save file
    Image_nuclei <- Image(image_nuclei)
    colorMode(Image_nuclei) <- "color"

    Image_nuclei <- paintObjects(x = nmask_watershed,
                                 tgt = Image_nuclei,
                                 col='#ff00ff')


    Image_nuclei_numbers <- Image(image_nuclei_numbers)
    colorMode(Image_nuclei_numbers) <- "color"

    Image_nuclei_numbers <- paintObjects(x = nmask_watershed,
                                        tgt = Image_nuclei_numbers,
                                        col='#ff00ff')

    # Display the number of nuclei
    print(paste("Number of nuclei: ", nucNo, sep=""))


    # Use the nucleus mask to leave only the nuclei part of the images
    Image_loaded <- EBImage::Image(data = image_loaded)
    Image_nucleus_part <- Image_loaded * EBImage::toRGB(nucleus_mask)

    # Use the nucleus mask to cut out nuclei of the images
    non_nucleus_mask <- 1-nucleus_mask
    Image_non_nucleus_part <- Image_loaded * EBImage::toRGB(non_nucleus_mask)

    # -------------------------------------------------------------------- #
    # -------------------- Statistics and data output -------------------- #
    # -------------------------------------------------------------------- #

    # Save results in data frame -------------------------------------------

    df_results[i,"dimension_x"] <- dim(image_loaded)[1]
    df_results[i,"dimension_y"] <- dim(image_loaded)[2]

    df_results[i,"intensity_sum_red_full"] <- sum(image_loaded[,,1])
    df_results[i,"intensity_sum_green_full"] <- sum(image_loaded[,,2])
    df_results[i,"intensity_sum_blue_full"] <- sum(image_loaded[,,3])

    df_results[i,"intensity_sum_red_nucleus_region"] <- sum(Image_nucleus_part[,,1])
    df_results[i,"intensity_sum_green_nucleus_region"] <- sum(Image_nucleus_part[,,2])
    df_results[i,"intensity_sum_blue_nucleus_region"] <- sum(Image_nucleus_part[,,3])

    df_results[i,"intensity_sum_red_without_nucleus_region"] <- sum(Image_non_nucleus_part[,,1])
    df_results[i,"intensity_sum_green_without_nucleus_region"] <- sum(Image_non_nucleus_part[,,2])
    df_results[i,"intensity_sum_blue_without_nucleus_region"] <- sum(Image_non_nucleus_part[,,3])

    if(!is.null(df_results)){
      write.csv(df_results,
                file = paste(output_dir, "intensity_summary.csv", sep=""), row.names = FALSE)

      write.csv2(df_results,
                 file = paste(output_dir, "intensity_summary_de.csv", sep=""), row.names = FALSE)
    }

    # Save all images ------------------------------------------------------

    # Normalized and histogram-adapted images
    tiff::writeTIFF(what = image_normalized,
                    where = paste(output_dir,
                                  image_name_wo_tif,
                                  "_normalized.tif",
                                  sep = ""),
                    bits.per.sample = 8L, compression = "none",
                    reduce = TRUE)

    tiff::writeTIFF(what = image_histogram_equalization,
                    where = paste(output_dir,
                                  image_name_wo_tif,
                                  "_histogram_equalized.tif",
                                  sep = ""),
                    bits.per.sample = 8L, compression = "none",
                    reduce = TRUE)

    # Images with marked nuclei
    tiff::writeTIFF(what = Image_nuclei,
                    where = paste(output_dir,
                                  image_name_wo_tif,
                                  "_nuclei.tif",
                                  sep = ""),
                    bits.per.sample = 8L, compression = "none",
                    reduce = TRUE)

    tiff::writeTIFF(what = Image_nuclei_numbers,
                    where = paste(output_dir,
                                  image_name_wo_tif,
                                  "_nuclei_numbers.tif",
                                  sep = ""),
                    bits.per.sample = 8L, compression = "none",
                    reduce = TRUE)

    # Images with left out nuclei and only with positions of nuclei
    tiff::writeTIFF(what = Image_nucleus_part,
                    where = paste(output_dir,
                                  image_name_wo_tif,
                                  "_nucleus_part.tif",
                                  sep = ""),
                    bits.per.sample = 8L, compression = "none",
                    reduce = TRUE)

    tiff::writeTIFF(what = Image_non_nucleus_part,
                    where = paste(output_dir,
                                  image_name_wo_tif,
                                  "_nucleus_left_out.tif",
                                  sep = ""),
                    bits.per.sample = 8L, compression = "none",
                    reduce = TRUE)

  } # end of the routine for every image in the directory

  return(df_results)

}

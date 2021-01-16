#' @title cellPixels
#' @description Main function for counting pixels in different regions
#' @details Input should be tif-format.
#' We are counting the intensities of pixels of different color within and
#' outside of the nuclei. We also add information about the entire image.
#' @aliases cellpixels CellPixels
#' @author Kai Budde
#' @export cellPixels
#' @param input_dir A character (directory that contains all images)
#' @param nucleus_color A character (color (layer) of nuclei)
#' @param protein_in_nuc_color A character (color (layer) of protein
#' expected in nucleus)
#' @param protein_in_cytosol_color A character (color (layer) of protein
#' expected in cytosol)
#' @param number_size_factor A number (factor to resize numbers for
#' numbering nuclei)
#' @param bit_depth A number (bit depth of the original czi image)
#' @param number_of_pixels_at_border_to_disregard A number (number of pixels
#' at the border of the image (rows and columns) that define the region
#' where found cells are disregarded)
#' @param add_scale_bar A logic (add scale bar to all images that are saved
#' if true)

cellPixels <- function(input_dir = NULL,
                       nucleus_color = "blue",
                       protein_in_nuc_color = "none",
                       protein_in_cytosol_color = "none",
                       number_size_factor = 0.2,
                       bit_depth = NULL,
                       number_of_pixels_at_border_to_disregard = 3,
                       add_scale_bar = FALSE) {

  # Basics and sourcing functions ------------------------------------------
  .old.options <- options()
  on.exit(options(.old.options))

  options(stringsAsFactors = FALSE, warn=-1)

  if(!("EBImage" %in% utils::installed.packages())){
    print("Installing EBImage.")
    BiocManager::install("EBImage")
  }

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
  #file_names <- file_names[grepl("tif$", file_names)]
  file_names <- file_names[grepl("czi$", file_names)]
  number_of_czis <- length(file_names)

  # Read in Python package for reading czi files
  # (Users will be asked to install miniconda
  # when starting for the first time)
  reticulate::py_install("czifile")
  zis <- reticulate::import("czifile")

  # If there is now tif-file, close function call
  if(number_of_czis == 0){
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
    "manual_quality_check" = rep(NA, number_of_czis),
    "dimension_x" = rep(NA, number_of_czis),
    "dimension_y" = rep(NA, number_of_czis),
    "number_of_nuclei" = rep(NA, number_of_czis),
    "color_of_second_protein_in_nuclei" = rep(NA, number_of_czis),
    "number_of_nuclei_with_second_protein" = rep(NA, number_of_czis),
    "color_of_third_protein_in_cytosol" = rep(NA, number_of_czis),
    "number_of_cells_with_third_protein" = rep(NA, number_of_czis),
    "intensity_sum_red_full" = rep(NA, number_of_czis),
    "intensity_sum_green_full" = rep(NA, number_of_czis),
    "intensity_sum_blue_full" = rep(NA, number_of_czis),
    "intensity_sum_red_nucleus_region" = rep(NA, number_of_czis),
    "intensity_sum_green_nucleus_region" = rep(NA, number_of_czis),
    "intensity_sum_blue_nucleus_region" = rep(NA, number_of_czis),
    "intensity_sum_red_without_nucleus_region" = rep(NA, number_of_czis),
    "intensity_sum_green_without_nucleus_region" = rep(NA, number_of_czis),
    "intensity_sum_blue_without_nucleus_region" = rep(NA, number_of_czis),
    "exposure_time_channel0" = rep(NA, number_of_czis),
    "exposure_time_channel1" = rep(NA, number_of_czis),
    "exposure_time_channel2" = rep(NA, number_of_czis))

  # Reduce the number of pixels for the borders because we will go from
  # 0 to number_of_pixels_at_border_to_disregard-1
  number_of_pixels_at_border_to_disregard <-
    number_of_pixels_at_border_to_disregard - 1


  # Go through every image in the directory --------------------------------
  for(i in 1:number_of_czis){

    print(paste("Dealing with file >>", file_names[i], "<<. (It is now ",
                Sys.time(), ".)", sep=""))

    ## Get the image name without the ".tif" ending
    #image_name_wo_czi <- gsub("\\.tif", "", file_names[i])
    # Get the image name without the ".czi" ending
    image_name_wo_czi <- gsub("\\.czi", "", file_names[i])

    # Get the image path
    image_path <- paste(input_dir, file_names[i], sep="/")

    ## Load image (hyperstack)
    #image_loaded <- tiff::readTIFF(source = image_path, info = FALSE,
    #                               all = TRUE,
    #                               convert = FALSE, as.is = TRUE)

    # Load image directly from czi and save bit depth of the image
    image_loaded <- zis$imread(image_path)

    czi_class <- zis$CziFile(image_path)
    bit_depth <- czi_class$dtype
    if(bit_depth == "uint16"){
      bit_depth <- 16
    }else if(bit_depth == "uint8"){
      bit_depth <- 8
    }else{
      print(paste("Something went wrong with the bit depth.", sep=""))
      return()
    }

    # ComponentBitCount -> shows the number of bits the camera can record
    # (is 0 or not specified if it is the same as dtype)
    metadata <- czi_class$metadata(czi_class)

    # Save the dimension information
    # B: (Acquisition) Block index in segmented experiments.
    # C: Channel in a Multi-Channel data set
    # X: Pixel index / offset in the X direction. Used for tiled images.
    # Y: Pixel index / offset in the y direction. Used for tiled images.
    # 0: Data contains uncompressed pixels.

    axes <- czi_class$axes
    axes <- unlist(strsplit(x = axes, split = ""))
    pos_channels <- grep(pattern = "C", x = axes)
    pos_x <- grep(pattern = "X", x = axes)
    pos_y <- grep(pattern = "Y", x = axes)

    camera_bit_depth <- gsub(
      pattern = ".+<ComponentBitCount>(.+)</ComponentBitCount>.+",
      replacement = "\\1",
      x = metadata)
    camera_bit_depth <- as.numeric(camera_bit_depth)
    if(!is.na(camera_bit_depth) && (camera_bit_depth != 0)){
      bit_depth <- camera_bit_depth
    }

    # Find red, green, and blue channel ID
    wavelengths <- gsub(
      pattern =  paste(".+<EmissionWavelength>(.+)</EmissionWavelength>.+",
                       ".+<EmissionWavelength>(.+)</EmissionWavelength>.+",
                       ".+<EmissionWavelength>(.+)</EmissionWavelength>.+",
                       sep=""),
      replacement = "\\1,\\2,\\3",
      x = metadata)
    wavelengths <- as.numeric(strsplit(wavelengths, split = ",")[[1]])

    red_id <- which(wavelengths == max(wavelengths))
    blue_id <- which(wavelengths == min(wavelengths))

    if(length(wavelengths) == 3){
      green_id <- c(1:3)[!(c(1:3) %in% red_id | c(1:3) %in% blue_id)]
    }else{
      print("There are more or fewer than 3 wavelengths.")
      return()
    }

    rgb_layers <- c(red_id, green_id, blue_id)

    rm(czi_class)

    # # Combine every channel of the image (separate list item) into one
    # # matrix
    # if(is.list(image_loaded)){
    #
    #   if(length(image_loaded) == 3){
    #     test_image <- array(dim = c(dim(image_loaded[[1]])[1],
    #                                 dim(image_loaded[[1]])[2],
    #                                 3))
    #     test_image[,,1] <- image_loaded[[1]]
    #     test_image[,,2] <- image_loaded[[2]]
    #     test_image[,,3] <- image_loaded[[3]]
    #
    #     image_loaded <- test_image
    #     rm(test_image)
    #   }else{
    #     print(paste("We do not have a tif-image with three layers",
    #                 " representing red, green, and blue.", sep=""))
    #     return()
    #   }
    #
    # }

    # Make a 3-D array of the image
    number_of_channels <- dim(image_loaded)[pos_channels]
    dim_x <- dim(image_loaded)[pos_x]
    dim_y <- dim(image_loaded)[pos_y]

    if(number_of_channels == 3){

      # Delete the dimensions of an array which have only one level
      image_loaded <- drop(image_loaded)

      # Permute dimensions of array
      pos_x <- pos_x - min(pos_x, pos_y, pos_channels) + 1
      pos_y <- pos_y - min(pos_x, pos_y, pos_channels) + 1
      pos_channels <- pos_channels - min(pos_x, pos_y, pos_channels) + 1

      image_loaded <- aperm(a = image_loaded, c(pos_y, pos_x, pos_channels))

      # Reorder the layers accoring to the colors
      copy_image_loaded <- image_loaded

      image_loaded[,,1] <- copy_image_loaded[,,rgb_layers[1]]
      image_loaded[,,2] <- copy_image_loaded[,,rgb_layers[2]]
      image_loaded[,,3] <- copy_image_loaded[,,rgb_layers[3]]

      rm(copy_image_loaded)

    }else{
      print(paste("We do not have a tif-image with three layers",
                  " representing red, green, and blue.", sep=""))
      return()
    }

    # Convert to intensities between 0 and 1

    # # Find the possible bit depth
    # if(is.null(bit_depth)){
    #   max_intensity <- max(image_loaded)
    #   if(max_intensity <= (2^8)-1){
    #     bit_depth <- 8
    #   }else if(max_intensity <= (2^12)-1){
    #     bit_depth <- 12
    #   }else if(max_intensity <= (2^16)-1){
    #     bit_depth <- 16
    #   }else{
    #     print(paste("Bit depth is higher than 16. Somehting might",
    #                 " be wrong.", sep=""))
    #     return()
    #   }
    #
    # }


    image_loaded <- image_loaded/(2^bit_depth-1)


    # Get the exposure time for every layer
    exposure_time <- gsub(
      pattern =  paste(".+<ExposureTime>(.+)</ExposureTime>.+",
                       ".+<ExposureTime>(.+)</ExposureTime>.+",
                       ".+<ExposureTime>(.+)</ExposureTime>.+",
                       sep=""),
      replacement = "\\1,\\2,\\3",
      x = metadata)

    exposure_time <- unlist(strsplit(x = exposure_time, split = ","))



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
    Image_nuclei <- EBImage::Image(image_nuclei)
    rm(image_nuclei)
    #display(Image_nuclei)

    # Blur the image
    Image_nuclei <- EBImage::gblur(Image_nuclei, sigma = 7)

    # Mask the nuclei
    if(grepl(pattern = "_20x_", file_names[i])){
      # Smaller moving rectangle if the magnification is 20x (instead of 40x)
      nmask <- EBImage::thresh(Image_nuclei, w=15, h=15, offset=0.01)
    }else{
      nmask <- EBImage::thresh(Image_nuclei, w=35, h=35, offset=0.01)
    }


    # Morphological opening to remove objects smaller than the structuring element
    nmask <- EBImage::opening(nmask, EBImage::makeBrush(5, shape='disc'))
    # Fill holes
    nmask <- EBImage::fillHull(nmask)

    # Save black-and-white-image as nucleus mask
    nucleus_mask <- nmask

    # Label each connected set of pixels with a distinct ID
    nmask <- EBImage::bwlabel(nmask)

    #display(nmask)

    # Record all nuclei that are at the edges of the image
    left  <- table(nmask[1:number_of_pixels_at_border_to_disregard,
                         1:dim(nmask)[2]])
    top   <- table(nmask[1:dim(nmask)[1],
                         1:number_of_pixels_at_border_to_disregard])
    right <- table(nmask[
      (dim(nmask)[1]-number_of_pixels_at_border_to_disregard):dim(nmask)[1],
      1:dim(nmask)[2]])
    bottom <- table(nmask[
      1:dim(nmask)[1],
      (dim(nmask)[2]-number_of_pixels_at_border_to_disregard):dim(nmask)[2]])


    left <- as.integer(names(left))
    top <- as.integer(names(top))
    right <- as.integer(names(right))
    bottom <- as.integer(names(bottom))

    nuclei_at_borders <- unique(c(left, top, right, bottom))

    # delete the 0 in the list (these are pixels that do not belong to a nucleus)
    nuclei_at_borders <- nuclei_at_borders[nuclei_at_borders != 0]

    # Delete all nuclei at border
    if(length(nuclei_at_borders) > 0){
      for(j in 1:length(nuclei_at_borders)){
        EBImage::imageData(nmask)[
          EBImage::imageData(nmask) == nuclei_at_borders[j]] <- 0
      }
      rm(j)
    }

    #display(nmask)

    # Delete all remaining nuclei that are smaller than 5% of the median size
    # object sizes
    # barplot(table(nmask)[-1])

    table_nmask <- table(nmask)
    nuc_min_size <- 0.1*stats::median(table_nmask[-1])

    # remove objects that are smaller than min_nuc_size
    to_be_removed <- as.integer(names(which(table_nmask < nuc_min_size)))
    if(length(to_be_removed) > 0){
      for(j in 1:length(to_be_removed)){
        EBImage::imageData(nmask)[
          EBImage::imageData(nmask) == to_be_removed[j]] <- 0
      }
      rm(j)
    }

    # Recount nuclei
    nmask <- EBImage::bwlabel(nmask)
    #display(nmask)

    # Watershed in order to distinct nuclei that are too close to each other
    nmask_watershed <-  EBImage::watershed(
      EBImage::distmap(nmask), tolerance = 0.5, ext = 1)
    #display(colorLabels(nmask_watershed), all=TRUE)

    # Count number of cells
    nucNo <- max(nmask_watershed)

    # Include numbers of nuclei
    Image_nuclei <- image_loaded
    Image_nuclei_numbers <- image_loaded

    table_nmask_watershed <- table(nmask_watershed)

    if(length(table_nmask_watershed[-1]) > 0){
      # remove 0
      nuc_numbers <- as.integer(names(table_nmask_watershed[-1]))

      for(j in 1:length(nuc_numbers)){

        # Find approximate midpoint of every nucleus
        dummy_coordinates <- which(
          EBImage::imageData(nmask_watershed) ==
            nuc_numbers[j], arr.ind = TRUE)


        pos_x <- round(mean(dummy_coordinates[,1]))
        pos_y <- round(mean(dummy_coordinates[,2]))

        Image_nuclei_numbers <- addNumberToImage(
          image = Image_nuclei_numbers,
          number = j,
          pos_x = pos_x,
          pos_y = pos_y,
          number_size_factor = number_size_factor,
          number_color = "green")
        Image_nuclei_numbers <- addNumberToImage(
          image = Image_nuclei_numbers,
          number = j,
          pos_x = pos_x,
          pos_y = pos_y,
          number_size_factor = number_size_factor,
          number_color = "blue")
      }
      rm(j)
    }

    # Add border of nuclei and save file
    Image_nuclei <- EBImage::Image(Image_nuclei)
    EBImage::colorMode(Image_nuclei) <- "color"

    Image_nuclei <- EBImage::paintObjects(
      x = nmask_watershed,
      tgt = Image_nuclei,
      thick = TRUE,
      col='#ff00ff')

    Image_nuclei_numbers <- EBImage::Image(Image_nuclei_numbers)
    EBImage::colorMode(Image_nuclei_numbers) <- "color"

    Image_nuclei_numbers <- EBImage::paintObjects(
      x = nmask_watershed,
      tgt = Image_nuclei_numbers,
      thick = TRUE,
      col='#ff00ff')

    # Display the number of nuclei
    print(paste("Number of nuclei: ", nucNo, sep=""))


    # Use the nucleus mask to leave only the nuclei part of the images
    Image_loaded <- EBImage::Image(data = image_loaded)
    Image_nucleus_part <- Image_loaded * EBImage::toRGB(nucleus_mask)

    # Use the nucleus mask to cut out nuclei of the images
    non_nucleus_mask <- 1-nucleus_mask
    Image_non_nucleus_part <- Image_loaded * EBImage::toRGB(non_nucleus_mask)


    # Count the number of nuclei that contain a second colored protein  ----

    if(!is.null(protein_in_nuc_color) & protein_in_nuc_color != "none"){

      # Save only color layer of the second protein colored
      image_protein_in_nuc <- getLayer(image = image_loaded,
                                       layer = protein_in_nuc_color)
      Image_protein_in_nuc <- EBImage::Image(image_protein_in_nuc)
      rm(image_protein_in_nuc)
      #display(Image_protein_in_nuc)

      # Blur the image
      Image_protein_in_nuc <- EBImage::gblur(Image_protein_in_nuc, sigma = 7)
      #display(Image_protein_in_nuc)


      # Mask the proteins within the nucleus
      if(grepl(pattern = "_20x_", file_names[i])){
        # Smaller moving rectangle if the magnification is 20x (instead of 40x)
        pmask <- EBImage::thresh(Image_protein_in_nuc, w=15, h=15, offset=0.07)
      }else{
        pmask <- EBImage::thresh(Image_protein_in_nuc, w=35, h=35, offset=0.07)
      }

      # Morphological opening to remove objects smaller than the structuring element
      pmask <- EBImage::opening(pmask, EBImage::makeBrush(5, shape='disc'))
      # Fill holes
      pmask <- EBImage::fillHull(pmask)
      #display(pmask)

      # Combine nmask and pmask and count the resulting nuclei containing
      # the proteins we are looking for

      n_p_mask <- nmask*pmask
      #display(n_p_mask)

      # Count the nuclei containing the proteins
      n_p_mask <- EBImage::bwlabel(n_p_mask)
      #display(n_p_mask)

      # Count number of cells
      nuc_with_proteins_No <- max(n_p_mask)

      # Add border of nuclei with proteins and save file
      Image_nuclei_numbers_proteins <- Image_nuclei_numbers
      EBImage::colorMode(Image_nuclei_numbers_proteins) <- "color"

      Image_nuclei_numbers_proteins <- EBImage::paintObjects(
        x = n_p_mask,
        tgt = Image_nuclei_numbers_proteins,
        opac = c(1,0.5),
        col=c('#FE1F14','#FE1F14'))

      # Display the number of nuclei with proteins
      print(paste("Number of nuclei that contain other colored proteins: ",
                  nuc_with_proteins_No, sep=""))

    }

    # Count the number of cells that contain a third colored protein  ------

    if(!is.null(protein_in_cytosol_color) & protein_in_cytosol_color != "none"){


      # Save only color layer of the second protein colored
      image_protein_in_cytosol <- getLayer(
        image = image_loaded, layer = protein_in_cytosol_color)

      Image_protein_in_cytosol <- EBImage::Image(image_protein_in_cytosol)
      rm(image_protein_in_cytosol)
      #display(Image_protein_in_cytosol)

      # Blur the image
      Image_protein_in_cytosol <- EBImage::gblur(Image_protein_in_cytosol, sigma = 5)
      #display(Image_protein_in_nuc)

      # Mask the proteins within the cytosol
      if(grepl(pattern = "_20x_", file_names[i])){
        # Smaller moving rectangle if the magnification is 20x (instead of 40x)
        cytosolmask <- EBImage::thresh(Image_protein_in_cytosol, w=100, h=100, offset=0.01)
      }else{
        cytosolmask <- EBImage::thresh(Image_protein_in_cytosol, w=200, h=200, offset=0.01)
      }

      #display(cytosolmask)

      # Morphological opening to remove objects smaller than the structuring element
      cytosolmask <- EBImage::opening(cytosolmask, EBImage::makeBrush(5, shape='disc'))
      # Fill holes
      cytosolmask <- EBImage::fillHull(cytosolmask)
      #display(cytosolmask)

      # Keep only those cell bodies that contain a nucleus
      cytosolmask <- EBImage::bwlabel(cytosolmask)
      no_of_cytosols <- max(cytosolmask)

      # Go through every stained cytosol and check for nucleus
      for(j in 1:no_of_cytosols){
        collocation_found <- max((cytosolmask==j)*nmask)
        if(collocation_found == 0){
          cytosolmask[cytosolmask==j] <- 0
        }
        rm(j)
      }

      cytosolmask <- EBImage::bwlabel(cytosolmask)
      no_of_cytosols <- max(cytosolmask)

      # Combine nmask and cytosolmask and count the resulting cell bodies
      # (nuclei that are within the stained proteins in the cytosol)
      n_c_mask <- nmask_watershed*(!cytosolmask==0)
      #display(n_c_mask)

      table_n_c_mask <- table(n_c_mask)

      # Count number of cells containing staining for protein
      cell_with_proteins_No <- length(names(table_n_c_mask))-1
      #display(n_c_mask)


      # remove objects that are smaller than min_nuc_size
      to_be_removed <- as.integer(names(which(table_n_c_mask < nuc_min_size)))
      if(length(to_be_removed) > 0){
        for(j in 1:length(to_be_removed)){
          EBImage::imageData(n_c_mask)[
            EBImage::imageData(n_c_mask) == to_be_removed[j]] <- 0
        }
        rm(j)
      }

      # Add border of cytosols with proteins and save file
      Image_cytosol_numbers_proteins <- Image_nuclei_numbers
      EBImage::colorMode(Image_cytosol_numbers_proteins) <- "color"

      Image_cytosol_numbers_proteins <- EBImage::paintObjects(
        x = cytosolmask,
        tgt = Image_cytosol_numbers_proteins,
        opac = c(1),
        col=c('#FFFF00'))

      #display(Image_cytosol_numbers_proteins)

      # Display the number of cells with proteins
      print(paste("Number of cells that contain other colored proteins: ",
                  cell_with_proteins_No, sep=""))

    }



    # -------------------------------------------------------------------- #
    # -------------------- Statistics and data output -------------------- #
    # -------------------------------------------------------------------- #

    # Save results in data frame -------------------------------------------

    df_results[i,"dimension_x"] <- dim(image_loaded)[2]
    df_results[i,"dimension_y"] <- dim(image_loaded)[1]
    df_results[i,"number_of_nuclei"] <- nucNo
    if(!is.null(protein_in_nuc_color) & protein_in_nuc_color != "none"){
      df_results[i,"color_of_second_protein_in_nuclei"] <- protein_in_nuc_color
      df_results[i,"number_of_nuclei_with_second_protein"] <- nuc_with_proteins_No
    }

    if(!is.null(protein_in_cytosol_color) & protein_in_cytosol_color != "none"){
      "color_of_third_protein_in_cytosol" <- protein_in_cytosol_color
      "number_of_cells_with_third_protein" <- cell_with_proteins_No
    }

    df_results[i,"intensity_sum_red_full"] <- sum(image_loaded[,,1])
    df_results[i,"intensity_sum_green_full"] <- sum(image_loaded[,,2])
    df_results[i,"intensity_sum_blue_full"] <- sum(image_loaded[,,3])

    df_results[i,"intensity_sum_red_nucleus_region"] <- sum(Image_nucleus_part[,,1])
    df_results[i,"intensity_sum_green_nucleus_region"] <- sum(Image_nucleus_part[,,2])
    df_results[i,"intensity_sum_blue_nucleus_region"] <- sum(Image_nucleus_part[,,3])

    df_results[i,"intensity_sum_red_without_nucleus_region"] <- sum(Image_non_nucleus_part[,,1])
    df_results[i,"intensity_sum_green_without_nucleus_region"] <- sum(Image_non_nucleus_part[,,2])
    df_results[i,"intensity_sum_blue_without_nucleus_region"] <- sum(Image_non_nucleus_part[,,3])

    df_results[i,"exposure_time_channel0"] <- exposure_time[1]
    df_results[i,"exposure_time_channel1"] <- exposure_time[2]
    df_results[i,"exposure_time_channel2"] <- exposure_time[3]

    # Save all images ------------------------------------------------------

    # Get information for adding a scale bar
    length_per_pixel_x <- gsub(
      pattern = paste(".+<Items>[[:space:]]+<Distance Id=\"X\">[[:space:]]+",
      "<Value>(.{2,16})</Value>.+",sep=""),
      replacement = "\\1",
      x = metadata)
    length_per_pixel_x <- tolower(length_per_pixel_x)
    length_per_pixel_x <- as.numeric(length_per_pixel_x)

    length_per_pixel_y <- gsub(
      pattern = paste(".+<Items>.+<Distance Id=\"Y\">[[:space:]]+",
                      "<Value>(.{2,16})</Value>.+",sep=""),
      replacement = "\\1",
      x = metadata)
    length_per_pixel_y <- tolower(length_per_pixel_y)
    length_per_pixel_y <- as.numeric(length_per_pixel_y)

    if(length_per_pixel_x != length_per_pixel_y){
      print("Dimension in x- and y-directions are different! ERROR!")
    }

    # Save metadata in txt file
    utils::write.table(metadata, file = paste(output_dir,image_name_wo_czi,
                                  "_metadata.txt", sep = ""),
                sep = "", row.names = FALSE, col.names = FALSE, quote = FALSE)


    # Original image (converted to tif)
    if(add_scale_bar){
      image_loaded <- addScaleBar(image = image_loaded,
                           length_per_pixel = length_per_pixel_x)
    }
    tiff::writeTIFF(what = image_loaded,
                    where = paste(output_dir,
                                  image_name_wo_czi,
                                  "_original.tif",
                                  sep = ""),
                    bits.per.sample = 8L, compression = "none",
                    reduce = TRUE)

    # Normalized and histogram-adapted images
    if(add_scale_bar){
      image_normalized <- addScaleBar(image = image_normalized,
                                  length_per_pixel = length_per_pixel_x)
      image_histogram_equalization <- addScaleBar(image = image_histogram_equalization,
                                                  length_per_pixel = length_per_pixel_x)
    }
    tiff::writeTIFF(what = image_normalized,
                    where = paste(output_dir,
                                  image_name_wo_czi,
                                  "_normalized.tif",
                                  sep = ""),
                    bits.per.sample = 8L, compression = "none",
                    reduce = TRUE)

    tiff::writeTIFF(what = image_histogram_equalization,
                    where = paste(output_dir,
                                  image_name_wo_czi,
                                  "_histogram_equalized.tif",
                                  sep = ""),
                    bits.per.sample = 8L, compression = "none",
                    reduce = TRUE)

    # Images with marked nuclei
    if(add_scale_bar){
      Image_nuclei <- addScaleBar(image = Image_nuclei,
                                  length_per_pixel = length_per_pixel_x)
      Image_nuclei_numbers <- addScaleBar(image = Image_nuclei_numbers,
                                          length_per_pixel = length_per_pixel_x)
    }
    tiff::writeTIFF(what = Image_nuclei,
                    where = paste(output_dir,
                                  image_name_wo_czi,
                                  "_nuclei.tif",
                                  sep = ""),
                    bits.per.sample = 8L, compression = "none",
                    reduce = TRUE)

    tiff::writeTIFF(what = Image_nuclei_numbers,
                    where = paste(output_dir,
                                  image_name_wo_czi,
                                  "_nuclei_numbers.tif",
                                  sep = ""),
                    bits.per.sample = 8L, compression = "none",
                    reduce = TRUE)

    # Images with marked nuclei and borders around the second protein in nuc
    if(add_scale_bar){
      Image_nuclei_numbers_proteins <- addScaleBar(
        image = Image_nuclei_numbers_proteins,
        length_per_pixel = length_per_pixel_x)
    }
    if(!is.null(protein_in_nuc_color) & protein_in_nuc_color != "none"){
      tiff::writeTIFF(what = Image_nuclei_numbers_proteins,
                      where = paste(output_dir,
                                    image_name_wo_czi,
                                    "_nuclei_numbers_proteins1_nuc.tif",
                                    sep = ""),
                      bits.per.sample = 8L, compression = "none",
                      reduce = TRUE)
    }

    # Images with marked nuclei and borders around the third protein in
    # cell bodies (proteins in cytosol)
    if(add_scale_bar){
      Image_cytosol_numbers_proteins <- addScaleBar(
        image = Image_cytosol_numbers_proteins,
        length_per_pixel = length_per_pixel_x)
    }
    if(!is.null(protein_in_cytosol_color) & protein_in_cytosol_color != "none"){
      tiff::writeTIFF(what = Image_cytosol_numbers_proteins,
                      where = paste(output_dir,
                                    image_name_wo_czi,
                                    "_nuclei_numbers_proteins2_cell.tif",
                                    sep = ""),
                      bits.per.sample = 8L, compression = "none",
                      reduce = TRUE)
    }

    # Images with left out nuclei and only with positions of nuclei
    if(add_scale_bar){
      Image_nucleus_part <- addScaleBar(
        image = Image_nucleus_part,
        length_per_pixel = length_per_pixel_x)
      Image_non_nucleus_part <- addScaleBar(
        image = Image_non_nucleus_part,
        length_per_pixel = length_per_pixel_x)
    }
    tiff::writeTIFF(what = Image_nucleus_part,
                    where = paste(output_dir,
                                  image_name_wo_czi,
                                  "_nucleus_part.tif",
                                  sep = ""),
                    bits.per.sample = 8L, compression = "none",
                    reduce = TRUE)

    tiff::writeTIFF(what = Image_non_nucleus_part,
                    where = paste(output_dir,
                                  image_name_wo_czi,
                                  "_nucleus_left_out.tif",
                                  sep = ""),
                    bits.per.sample = 8L, compression = "none",
                    reduce = TRUE)

    # Remove all variables but the ones used before the for loop

    list_of_variables <- ls()
    keep_variables <- c("df_results", "zis",
                        "bit_depth", "file_names", "input_dir",
                        "nucleus_color","number_of_czis",
                        "number_of_pixels_at_border_to_disregard",
                        "number_size_factor", "output_dir",
                        "protein_in_cytosol_color",
                        "protein_in_nuc_color",
                        "add_scale_bar",
                        ".old.options")

    remove_variables <- list_of_variables[
      !(list_of_variables %in% keep_variables)]

    rm(list = remove_variables)
    rm(remove_variables)


  } # end of the routine for every image in the directory

  if(!is.null(df_results)){
    utils::write.csv(df_results,
                     file = paste(output_dir, "intensity_summary.csv", sep=""), row.names = FALSE)

    utils::write.csv2(df_results,
                      file = paste(output_dir, "intensity_summary_de.csv", sep=""), row.names = FALSE)
  }

  return(df_results)

}

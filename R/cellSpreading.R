#' @title cellSpreading
#' @description Main function for counting pixels in different regions
#' @details Input should be czi or tif-format with dim(z)>=1.
#' We are trying to identify different cells by using nuclei and Actin layers.
#' @aliases cellspreading CellSpreading Cellspreading
#' @author Kai Budde-Sagert
#' @export cellSpreading
#' @param input_dir A character (directory that contains all images)
#' @param nucleus_color A character (color (layer) of nuclei)
#' @param cell_body_color A character (color (layer) of cell body)
#' @param number_of_pixels_at_border_to_disregard A number (number of pixels
#' at the border of the image (rows and columns) that define the region
#' where found cells are disregarded)
#' @param use_histogram_equalized A boolean (use histogram equalized images for everything)
#' @param normalize_nuclei_layer A boolean (state whether nucleus layer should be normalized)
#' @param magnification_objective A number (magnification of objective if not given in metadata or if metadata is wrong)
#' @param apotome_section A boolean (TRUE is sectioned image shall be used)
#' @param blur_sigma A number (blurring factor)
#' @param thresh_w_h_nuc A number (width and heigth for thresholding)
#' @param thresh_offset A number (offset for thresholding)
#' @param number_size_factor A number (factor to resize numbers for
#' numbering nuclei)
#' @param add_scale_bar A logic (add scale bar to all images that are saved
#' if true)

cellSpreading <- function(input_dir = NULL,
                          nucleus_color = "blue",
                          cell_body_color = "red",
                          number_of_pixels_at_border_to_disregard = 3,
                          use_histogram_equalized = FALSE,
                          normalize_nuclei_layer = FALSE,
                          magnification_objective = NULL,
                          apotome_section = FALSE,
                          blur_sigma = NULL,
                          thresh_w_h_nuc = NULL,
                          thresh_offset = NULL,
                          number_size_factor = 1,
                          add_scale_bar = FALSE){

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

  # Save the file names (prefer czi, if available) -------------------------
  file_names <- list.files(path = input_dir)
  file_names <- file_names[grepl("czi$", file_names)]

  if(length(file_names) == 0){
    print("No czi files found, looking for tif files.")
    file_names <- list.files(path = input_dir)
    file_names <- file_names[grepl("tif$", file_names)]
    image_format <- "tif"
  }else{
    image_format <- "czi"
  }

  number_of_images <- length(file_names)

  # # Read in Python package for reading czi files
  # # (Users will be asked to install miniconda
  # # when starting for the first time)
  # reticulate::py_install("czifile")
  # zis <- reticulate::import("czifile")

  # If there is now image file (czi or tif), close function call
  if(number_of_images == 0){
    print("No image found.")
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
    "manual_quality_check" = rep(NA, number_of_images),
    "dimension_x" = rep(NA, number_of_images),
    "dimension_y" = rep(NA, number_of_images),
    "number_of_nuclei" = rep(NA, number_of_images),
    "nuclei_area_mean_in_pixels" = rep(NA, number_of_images),
    "nuclei_area_sd_in_pixels" = rep(NA, number_of_images),
    "cell_area_mean_in_pixels" = rep(NA, number_of_images),
    "cell_area_sd_in_pixels" = rep(NA, number_of_images))

  # Reduce the number of pixels for the borders because we will go from
  # 0 to number_of_pixels_at_border_to_disregard-1
  number_of_pixels_at_border_to_disregard <-
    number_of_pixels_at_border_to_disregard - 1


  # Go through every image in the directory --------------------------------
  for(i in 1:number_of_images){

    print(paste("Dealing with file >>", file_names[i], "<<. (It is now ",
                Sys.time(), ".)", sep=""))

    if(image_format == "czi"){
      # Get the image name without the ".czi" ending
      image_name_wo_czi <- gsub("\\.czi", "", file_names[i])

    }else if(image_format == "tif"){
      # Get the image name without the ".tif" ending
      image_name_wo_czi <- gsub("\\.tif", "", file_names[i])

    }

    # Get the image path
    image_path <- paste(input_dir, file_names[i], sep="/")


    # Load image (and metainformation if available)

    if(image_format == "czi"){
      # Load image directly from czi (as EBImage)
      #image_loaded <- zis$imread(image_path)
      image_loaded <- readCzi::readCzi(input_file = image_path)

      # Read metadata
      df_metadata <- readCzi::readCziMetadata(input_file = image_path,
                                              save_metadata = FALSE)

      # Stack image if it contains a zstack
      if(df_metadata$dim_z > 1){
        image_loaded <- readCzi::stackLayers(image_data = image_loaded,
                                             stack_method = "average",
                                             as_array = TRUE)
      }else{

        if(dim(image_loaded)[4] != df_metadata$dim_z){
          print("Something went wrong with the z-dimension.")
          return()
        }

        # Drop z-layer of multidimensional array if dim_z == 1

        image_loaded <- drop(image_loaded)
      }

    }else if(image_format == "tif"){
      # Load image
      # image_loaded <- tiff::readTIFF(source = image_path, info = FALSE,
      #                                all = TRUE,
      #                                convert = FALSE, as.is = FALSE)
      image_loaded <- EBImage::readImage(files = image_path, type = "tiff")

      # Read metadata
      df_metadata <- read.csv(file = metadata_file)
      df_metadata <- df_metadata[
        gsub(pattern = "\\.czi", replacement = "", x = df_metadata$fileName) ==
          gsub(pattern = "_co\\.tif", replacement = "", x = file_names[i]),]

      # Reorder channels
      image_dummy <- image_loaded
      if(!is.na(df_metadata$red_channel[1])){
        image_loaded[,,1] <- image_dummy[,,df_metadata$red_channel[1]]
      }
      if(!is.na(df_metadata$red_channel[1])){
        image_loaded[,,2] <- image_dummy[,,df_metadata$green_channel[1]]
      }
      if(!is.na(df_metadata$red_channel[1])){
        image_loaded[,,3] <- image_dummy[,,df_metadata$blue_channel[1]]
      }

    }

    # Number of total pixels
    number_of_total_pixels <- dim(image_loaded)[1]*dim(image_loaded)[2]


    # -------------------------------------------------------------------- #
    # --------------------------- Find Nuclei ---------------------------- #
    # -------------------------------------------------------------------- #

    # Normalize intensity --------------------------------------------------
    image_normalized <- readCzi::normalizeIntensity(image = image_loaded)

    # Use Contrast Limited Adaptive Histogram Equalization
    image_histogram_equalization <- EBImage::clahe(x = image_loaded, nx = 4)

    # Save only color layer of nuclei
    if(use_histogram_equalized){
      image_nuclei <- getLayer(image = image_histogram_equalization, layer = nucleus_color)
    }else{
      image_nuclei <- getLayer(image = image_loaded, layer = nucleus_color)
    }

    # Brighten nuclei image if required
    if(normalize_nuclei_layer || apotome_section){
      image_nuclei <- image_nuclei/max(image_nuclei)
    }

    Image_nuclei <- EBImage::Image(image_nuclei)
    rm(image_nuclei)
    #display(Image_nuclei)

    if(is.null(magnification_objective)){
      magnification <- df_metadata$objective_magnification[df_metadata$fileName == file_names[i]]
    }else{
      magnification <- magnification_objective
    }

    # Blur the image
    if(is.null(blur_sigma) | apotome_section){
      if(apotome_section){
        Image_nuclei <- EBImage::gblur(Image_nuclei, sigma = 20) #14
      }else{
        # Image_nuclei <- EBImage::gblur(Image_nuclei, sigma = 7)
        Image_nuclei <- EBImage::gblur(Image_nuclei, sigma = 6)
      }
    }else if(blur_sigma > 0){
      Image_nuclei <- EBImage::gblur(Image_nuclei, sigma = blur_sigma)
    }

    # Mask the nuclei
    if(is.null(thresh_w_h_nuc) || is.null(thresh_offset)){

      if(magnification < 20){
        # Smaller moving rectangle if the objective magnification is e.g. 10x
        nmask <- EBImage::thresh(Image_nuclei, w=8, h=8, offset=0.01)
      }else if(magnification < 40){ # Magnification of objective could be 20
        # Smaller moving rectangle if the objective magnification is e.g. 20x
        nmask <- EBImage::thresh(Image_nuclei, w=15, h=15, offset=0.01)
      }else if(magnification < 60){
        # Objective magnification of 40x
        nmask <- EBImage::thresh(Image_nuclei, w=30, h=30, offset=0.01)
      }else{
        # Objective magnification of 63x
        if(use_histogram_equalized){
          nmask <- EBImage::thresh(Image_nuclei, w=dim(Image_nuclei)[1]/4, h=dim(Image_nuclei)[1]/4, offset=0.05)
        }else{
          nmask <- EBImage::thresh(Image_nuclei, w=dim(Image_nuclei)[1]/4, h=dim(Image_nuclei)[1]/4, offset=0.02)
        }
      }
    }else{
      nmask <- EBImage::thresh(Image_nuclei, w=thresh_w_h_nuc, h=thresh_w_h_nuc, offset=thresh_offset)
    }
    # display(nmask)


    # Morphological opening to remove objects smaller than the structuring element
    nmask <- EBImage::opening(nmask, EBImage::makeBrush(5, shape='disc'))
    # Fill holes
    nmask <- EBImage::fillHull(nmask)

    # Label each connected set of pixels with a distinct ID
    nmask <- EBImage::bwlabel(nmask)

    #display(nmask)

    if(magnification < 20){
      nmask_watershed <-  EBImage::watershed(
        EBImage::distmap(nmask), tolerance = 0.45, ext = 3) # it was tolerance = 0.5
    }else if(magnification < 40){
      nmask_watershed <-  EBImage::watershed(
        EBImage::distmap(nmask), tolerance = 0.1, ext = 6) #sqrt(15), see thresh(Image_nuclei...
    }else if(magnification < 60){
      nmask_watershed <-  EBImage::watershed(
        EBImage::distmap(nmask), tolerance = 0.1, ext = 6) # used to be 045, 9 | 0.3, 3
    }else{
      nmask_watershed <-  EBImage::watershed(
        EBImage::distmap(nmask), tolerance = 0.3, ext = 6) # used to be 045, 9 | 0.3, 3
    }

    #display(colorLabels(nmask_watershed), all=TRUE)
    #display(colorLabels(nmask), all=TRUE)
    # display(colorLabels(EBImage::watershed(EBImage::distmap(nmask), tolerance = 0.2, ext = 6)), all=TRUE)

    nmask <- nmask_watershed


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

    # display(colorLabels(nmask), all=TRUE)

    # Delete all remaining nuclei that are smaller than 5% of the median size
    # object sizes
    # barplot(table(nmask)[-1])

    table_nmask <- table(nmask)
    nuc_min_size <- 0.2*stats::median(table_nmask[-1])

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
    table_nmask <- table(nmask)
    number_of_nuclei <- length(table_nmask)-1
    if(number_of_nuclei >= 1){
      for(nuc_id in 1:number_of_nuclei){
        current_nuc_id <- as.numeric(names(table_nmask))[nuc_id+1]
        nmask[nmask == current_nuc_id] <- nuc_id
      }
    }

    # Count number of cells
    nmask_watershed <- nmask
    # display(colorLabels(nmask_watershed), all=TRUE)

    nucNo <- max(nmask_watershed)

    # Mean and standard deviation of nuclei masks
    df_nmask <- as.data.frame(table(nmask_watershed))
    names(df_nmask) <- c("ID", "Frequency")
    df_nmask <- df_nmask[!df_nmask$ID == 0,]
    nuc_mean <- mean(df_nmask$Frequency)
    nuc_sd <- sd(df_nmask$Frequency)

    # Include numbers of nuclei
    Image_nuclei <- image_loaded
    Image_nuclei_numbers <- as.array(image_loaded)

    table_nmask_watershed <- table(nmask_watershed)

    # Sort the nuclei from top to bottom and left to right
    if(length(table_nmask_watershed[-1]) > 0){
      df_nmask_watershed <- data.frame(table(nmask_watershed))

      names(df_nmask_watershed)[names(df_nmask_watershed) == "nmask_watershed"] <- "NucNo"
      df_nmask_watershed <- df_nmask_watershed[!df_nmask_watershed$NucNo == 0,]

      df_nmask_watershed$x_coord <- NA
      df_nmask_watershed$y_coord <- NA
      for(j in 1:nrow(df_nmask_watershed)){

        # Find approximate midpoint of every nucleus
        dummy_coordinates <- which(
          EBImage::imageData(nmask_watershed) ==
            j, arr.ind = TRUE)


        pos_x <- round(mean(dummy_coordinates[,1]))
        pos_y <- round(mean(dummy_coordinates[,2]))


        df_nmask_watershed$x_coord[j] <- pos_x
        df_nmask_watershed$y_coord[j] <- pos_y

      }
      rm(j)

      df_nmask_watershed <- df_nmask_watershed[order(df_nmask_watershed$y_coord, df_nmask_watershed$x_coord),]
      rownames(df_nmask_watershed) <- NULL
      df_nmask_watershed$newNucNo <- 1:nrow(df_nmask_watershed)

      for(j in 1:nrow(df_nmask_watershed)){

        Image_nuclei_numbers <- addNumberToImage(
          image = Image_nuclei_numbers,
          number = df_nmask_watershed$newNucNo[j],
          pos_x = df_nmask_watershed$x_coord[j],
          pos_y = df_nmask_watershed$y_coord[j],
          number_size_factor = number_size_factor,
          number_color = nucleus_color)
        Image_nuclei_numbers <- addNumberToImage(
          image = Image_nuclei_numbers,
          number = df_nmask_watershed$newNucNo[j],
          pos_x = df_nmask_watershed$x_coord[j],
          pos_y = df_nmask_watershed$y_coord[j],
          number_size_factor = number_size_factor,
          number_color = cell_body_color)
      }
      rm(j)

    }

    # Add border of nuclei and save file
    Image_nuclei <- EBImage::Image(Image_nuclei, colormode = "color")
    # EBImage::colorMode(Image_nuclei) <- "color"

    Image_nuclei <- EBImage::paintObjects(
      x = nmask_watershed,
      tgt = Image_nuclei,
      thick = TRUE,
      col='#ff00ff')

    Image_nuclei_numbers <- EBImage::Image(Image_nuclei_numbers, colormode = "color")
    # EBImage::colorMode(Image_nuclei_numbers) <- "color"

    Image_nuclei_numbers <- EBImage::paintObjects(
      x = nmask_watershed,
      tgt = Image_nuclei_numbers,
      thick = TRUE,
      col='#ff00ff')

    # Display the number of nuclei
    print(paste("Number of nuclei: ", nucNo, sep=""))


    # -------------------------------------------------------------------- #
    # ------------------------ Find cell bodies -------------------------- #
    # -------------------------------------------------------------------- #


    # Save only color layer of the cell bodies
    if(use_histogram_equalized){
      image_cell_body <- getLayer(image = image_histogram_equalization, layer = cell_body_color)
    }else{
      image_cell_body <- getLayer(image = image_loaded, layer = cell_body_color)
    }

    # Brighten cell body image if required
    if(normalize_nuclei_layer || apotome_section){
      image_cell_body <- image_cell_body/max(image_cell_body)
    }

    Image_cell_body <- EBImage::Image(image_cell_body)
    rm(image_cell_body)
    #display(Image_cell_body)

    # Blur the image
    if(is.null(blur_sigma) | apotome_section){
      if(apotome_section | magnification == 63){
        Image_cell_body <- EBImage::gblur(Image_cell_body, sigma = 20)
      }else{
        Image_cell_body <- EBImage::gblur(Image_cell_body, sigma = 6)
      }
    }else if(blur_sigma > 0){
      Image_cell_body <- EBImage::gblur(Image_cell_body, sigma = blur_sigma)
    }

    #display(Image_cell_body)

    # Mask the cell bodies
    if(is.null(thresh_w_h_nuc) || is.null(thresh_offset)){

      if(magnification < 20){
        # Smaller moving rectangle if the objective magnification is e.g. 10x
        cellmask <- EBImage::thresh(Image_cell_body, w=8, h=8, offset=0.01)
      }else if(magnification < 40){ # Magnification of objective could be 20
        # Smaller moving rectangle if the objective magnification is e.g. 20x
        cellmask <- EBImage::thresh(Image_cell_body, w=15, h=15, offset=0.01)
      }else if(magnification < 60){
        # Objective magnification of 40x
        cellmask <- EBImage::thresh(Image_cell_body, w=30, h=30, offset=0.01)
      }else{
        # Objective magnification of 63x
        cellmask <- EBImage::thresh(Image_cell_body, w=dim(Image_nuclei)[1]/4, h=dim(Image_nuclei)[1]/4, offset=0.01)
      }
    }else{
      cellmask <- EBImage::thresh(Image_cell_body, w=thresh_w_h_nuc, h=thresh_w_h_nuc, offset=thresh_offset)
    }
    # display(cellmask)



    # Morphological opening to remove objects smaller than the structuring element
    cellmask <- EBImage::opening(cellmask, EBImage::makeBrush(5, shape='disc'))
    # Fill holes
    cellmask <- EBImage::fillHull(cellmask)

    # Label each connected set of pixels with a distinct ID
    cellmask <- EBImage::bwlabel(cellmask)

    # Remove all areas without a nuclei mask in it
    cells <- table(cellmask)
    cells <- as.numeric(names(cells))
    cells <- cells[!cells == 0]


    for(j in cells){
      nucleus_in_cell <- (cellmask == j) * nmask_watershed
      display(nucleus_in_cell)

      if(sum(nucleus_in_cell) == 0){
        EBImage::imageData(cellmask)[cellmask == j] <- 0
      }
    }
    rm(j)

    # display(cellmask)

    # Segmentation of cell bodies
    cellmask <- EBImage::propagate(Image_cell_body, seeds=nmask, mask=cellmask, lambda = 1e-8)

    #display(colorLabels(cellmask), all=TRUE)

    # Mean and standard deviation of nuclei masks
    df_cellmask <- as.data.frame(table(cellmask))
    names(df_cellmask) <- c("ID", "Frequency")
    df_cellmask <- df_cellmask[!df_cellmask$ID == 0,]
    cell_mean <- mean(df_cellmask$Frequency)
    cell_sd <- sd(df_cellmask$Frequency)


    # Add border of cell bodies and save file
    Image_cellbodies <- EBImage::Image(Image_nuclei_numbers, colormode = "color")

    Image_cellbodies <- EBImage::paintObjects(
      x = cellmask,
      tgt = Image_cellbodies,
      thick = TRUE,
      col='#ffff00')

    # display(Image_cellbodies)

    # -------------------------------------------------------------------- #
    # -------------------- Statistics and data output -------------------- #
    # -------------------------------------------------------------------- #

    # Save results in data frame -------------------------------------------

    df_results[i,"dimension_x"] <- dim(image_loaded)[2]
    df_results[i,"dimension_y"] <- dim(image_loaded)[1]
    df_results[i,"number_of_nuclei"] <- nucNo
    df_results[i,"nuclei_area_mean_in_pixels"] <- nuc_mean
    df_results[i,"nuclei_area_sd_in_pixels"] <- nuc_sd
    df_results[i,"cell_area_mean_in_pixels"] <- cell_mean
    df_results[i,"cell_area_sd_in_pixels"] <- cell_sd

    # Save all images ------------------------------------------------------

    length_per_pixel_x_in_um <- df_metadata$scaling_x_in_um[1]
    length_per_pixel_y_in_um <- df_metadata$scaling_y_in_um[1]

    if(round(length_per_pixel_x_in_um, digits = 5) !=
       round(length_per_pixel_y_in_um, digits = 5)){
      print("Dimension in x- and y-directions are different! ERROR!")
    }

    # Original image (converted to tif)
    if(add_scale_bar){
      image_loaded <- addScaleBar(image = image_loaded,
                                  length_per_pixel = length_per_pixel_x_in_um)
    }

    Image_loaded <- EBImage::Image(data = image_loaded, colormode = "Color")
    EBImage::writeImage(x = Image_loaded,
                        files = paste(output_dir, image_name_wo_czi, "_original.tif", sep = ""),
                        type = "tiff",
                        bits.per.sample = 8)

    # Normalized and histogram-adapted images
    if(add_scale_bar){
      image_normalized <- addScaleBar(image = image_normalized,
                                      length_per_pixel = length_per_pixel_x_in_um)
      image_histogram_equalization <- addScaleBar(image = image_histogram_equalization,
                                                  length_per_pixel = length_per_pixel_x_in_um)
    }

    Image_normalized <- EBImage::Image(data = image_normalized, colormode = "Color")
    EBImage::writeImage(x = Image_normalized,
                        files = paste(output_dir, image_name_wo_czi, "_normalized.tif", sep = ""),
                        type = "tiff",
                        bits.per.sample = 8)

    Image_histogram_equalization <- EBImage::Image(data = image_histogram_equalization, colormode = "Color")
    EBImage::writeImage(x = Image_histogram_equalization,
                        files = paste(output_dir, image_name_wo_czi, "_histogram_equalized.tif", sep = ""),
                        type = "tiff",
                        bits.per.sample = 8)

    # Images with marked nuclei
    if(add_scale_bar){
      Image_nuclei <- addScaleBar(image = Image_nuclei,
                                  length_per_pixel = length_per_pixel_x_in_um)
      Image_nuclei_numbers <- addScaleBar(image = Image_nuclei_numbers,
                                          length_per_pixel = length_per_pixel_x_in_um)
    }

    EBImage::writeImage(x = Image_nuclei,
                        files = paste(output_dir, image_name_wo_czi, "_nuclei.tif", sep = ""),
                        type = "tiff",
                        bits.per.sample = 8)

    EBImage::writeImage(x = Image_nuclei_numbers,
                        files = paste(output_dir, image_name_wo_czi, "_nuclei_numbers.tif", sep = ""),
                        type = "tiff",
                        bits.per.sample = 8)


    # Images with marked nuclei and marked cell bodies
    if(add_scale_bar){
      Image_cellbodies <- addScaleBar(image = Image_cellbodies,
                                  length_per_pixel = length_per_pixel_x_in_um)
    }

    EBImage::writeImage(x = Image_cellbodies,
                        files = paste(output_dir, image_name_wo_czi, "_cellbodies.tif", sep = ""),
                        type = "tiff",
                        bits.per.sample = 8)


    # Remove all variables but the ones used before the for loop

    list_of_variables <- ls()
    keep_variables <- c("input_dir", "nucleus_color", "cell_body_color",
                        "number_of_pixels_at_border_to_disregard",
                        "use_histogram_equalized",
                        "normalize_nuclei_layer",
                        "magnification_objective",
                        "apotome_section", "blur_sigma",
                        "thresh_w_h_nuc", "thresh_offset",
                        "number_size_factor", "add_scale_bar",
                        "df_results", "file_names", "image_format",
                         "output_dir", ".old.options")


    remove_variables <- list_of_variables[
      !(list_of_variables %in% keep_variables)]

    rm(list = remove_variables)
    rm(remove_variables)


  } # end of the routine for every image in the directory
  rm(i)

  if(!is.null(df_results)){
    readr::write_csv(df_results,
                     file = paste(output_dir, "image_analysis_summary_en.csv", sep=""), row.names = FALSE)

    readr::write_csv2(df_results,
                      file = paste(output_dir, "image_analysis_summary_de.csv", sep=""), row.names = FALSE)
  }

  return(df_results)

}

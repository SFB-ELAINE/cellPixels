#' @title cellPixels
#' @description Main function for counting pixels in different regions
#' @details Input should be czi or tif-format with dim(z)>=1.
#' We are counting the intensities of pixels of different color within and
#' outside of the nuclei. We also add information about the entire image.
#' @aliases cellpixels CellPixels
#' @author Kai Budde
#' @export cellPixels
#' @param input_dir A character (directory that contains all images)
#' @param apotome A boolean (TRUE if Apotome was used)
#' @param apotome_section A boolean (TRUE is sectioned image shall be used)
#' @param nucleus_color A character (color (layer) of nuclei)
#' @param protein_in_nuc_color A character (color (layer) of protein
#' expected in nucleus)
#' @param protein_in_cytosol_color A character (color (layer) of protein
#' expected in cytosol)
#' @param protein_in_membrane_color A character (color (layer) of protein
#' expected in membrane)
#' @param number_size_factor A number (factor to resize numbers for
#' numbering nuclei)
#' @param bit_depth A number (bit depth of the original czi image)
#' @param number_of_pixels_at_border_to_disregard A number (number of pixels
#' at the border of the image (rows and columns) that define the region
#' where found cells are disregarded)
#' @param add_scale_bar A logic (add scale bar to all images that are saved
#' if true)
#' @param thresh_w_h_nuc A number (right now, both thresh_w_h_nuc and
#' thresh_offset must be supplied even if only one of them should be changed
#' manually)
#' @param thresh_offset A number (threshold for finding nuclei)
#' @param thresh_offset_protein_in_nucleus A number (threshold for
#' determining "positive" cells containing a certain protein in nucleus area)
#' @param thresh_offset_protein_in_cytosol A number (threshold for
#' determining "positive" cells containing a certain protein in cytosol area)
#' @param metadata_file A character (file with meta data if tifs are used)
#' @param normalize_nuclei_layer A boolean (state whether nucleus layer should be normalized)
#' @param magnification_objective A number (magnification of objective if not given in metadata or if metadata is wrong)

cellPixels <- function(input_dir = NULL,
                       apotome = FALSE,
                       apotome_section = FALSE,
                       nucleus_color = "blue",
                       protein_in_nuc_color = "none",
                       protein_in_cytosol_color = "none",
                       protein_in_membrane_color = "none",
                       number_size_factor = 0.2,
                       bit_depth = NULL,
                       number_of_pixels_at_border_to_disregard = 3,
                       add_scale_bar = FALSE,
                       thresh_w_h_nuc = NULL,
                       thresh_offset = NULL,
                       thresh_offset_protein_in_nucleus = NULL,
                       thresh_offset_protein_in_cytosol = NULL,
                       blur_sigma = NULL,
                       use_histogram_equalized = FALSE,
                       metadata_file = NULL,
                       normalize_nuclei_layer = FALSE,
                       magnification_objective = NULL) {

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
    "color_of_second_protein_in_nuclei" = rep(NA, number_of_images),
    "number_of_nuclei_with_second_protein" = rep(NA, number_of_images),
    "color_of_third_protein_in_cytosol" = rep(NA, number_of_images),
    "number_of_cells_with_third_protein" = rep(NA, number_of_images),
    "intensity_sum_red_full" = rep(NA, number_of_images),
    "intensity_sum_green_full" = rep(NA, number_of_images),
    "intensity_sum_blue_full" = rep(NA, number_of_images),
    "intensity_sum_red_nucleus_region" = rep(NA, number_of_images),
    "intensity_sum_green_nucleus_region" = rep(NA, number_of_images),
    "intensity_sum_blue_nucleus_region" = rep(NA, number_of_images),
    "intensity_sum_red_without_nucleus_region" = rep(NA, number_of_images),
    "intensity_sum_green_without_nucleus_region" = rep(NA, number_of_images),
    "intensity_sum_blue_without_nucleus_region" = rep(NA, number_of_images),
    "intensity_sum_red_foreground" = rep(NA, number_of_images),
    "intensity_sum_green_foreground" = rep(NA, number_of_images),
    "intensity_sum_blue_foreground" = rep(NA, number_of_images),
    "intensity_sum_red_foreground_without_nucleus_region" = rep(NA, number_of_images),
    "intensity_sum_green_foreground_without_nucleus_region" = rep(NA, number_of_images),
    "intensity_sum_blue_foreground_without_nucleus_region" = rep(NA, number_of_images),
    "intensity_mean_red_background" = rep(NA, number_of_images),
    "intensity_mean_green_background" = rep(NA, number_of_images),
    "intensity_mean_blue_background" = rep(NA, number_of_images),
    "number_of_total_pixels" = rep(NA, number_of_images),
    "number_of_pixels_nucleus_region" = rep(NA, number_of_images),
    "number_of_pixels_foreground" = rep(NA, number_of_images),
    "number_of_pixels_foreground_without_nucleus_region" = rep(NA, number_of_images),
    "number_of_clusters" = rep(NA, number_of_images),
    "mean_cluster_size" = rep(NA, number_of_images),
    "median_cluster_size" = rep(NA, number_of_images),
    "image_cropped" = rep(NA, number_of_images)
    )

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

      # czi_class <- zis$CziFile(image_path)
      # bit_depth <- czi_class$dtype
      # if(bit_depth == "uint16"){
      #   bit_depth <- 16
      # }else if(bit_depth == "uint8"){
      #   bit_depth <- 8
      # }else{
      #   print(paste("Something went wrong with the bit depth.", sep=""))
      #   return()
      # }

      # ComponentBitCount -> shows the number of bits the camera can record
      # (is 0 or not specified if it is the same as dtype)
      # metadata <- czi_class$metadata(czi_class)

      # Save the dimension information
      # X: Pixel index / offset in the X direction. Used for tiled images.
      # Y: Pixel index / offset in the y direction. Used for tiled images.
      # C: Channel in a Multi-Channel data set.
      # Z: Slice index (Z – direction).
      # T: Time point in a sequentially acquired series of data.
      # R: Rotation – used in acquisition modes where the data is recorded
      #    from various angles.
      # S: Scene – for clustering items in X/Y direction (data belonging to
      #   contiguous regions of interests in a mosaic image).
      # B: (Acquisition) Block index in segmented experiments.
      # H: Phase index – for specific acquisition methods.
      # M: Mosaic tile index – to reconstruct the  image acquisition order
      #    and to be used in conjunction with a global position list and the
      #    scaling information to define the pixel offset (alternative to
      #    using StartX, StartY to allow multi-M SubBlocks in case of
      #    homogenous mosaics, i.e. all tile share the same position list)

      # 0: Data contains uncompressed pixels.

      # axes <- czi_class$axes
      # axes <- unlist(strsplit(x = axes, split = ""))
      # pos_channels <- grep(pattern = "C", x = axes)
      # number_of_channels <- dim(image_loaded)[pos_channels]
      # pos_x <- grep(pattern = "X", x = axes)
      # dim_x <- dim(image_loaded)[pos_x]
      # pos_y <- grep(pattern = "Y", x = axes)
      # dim_y <- dim(image_loaded)[pos_y]
      #
      # # Check for multiple scenes
      # if("S" %in% axes){
      #   pos_scenes <- grep(pattern = "S", x = axes)
      #   number_scenes <- dim(image_loaded)[pos_scenes]
      #   if(number_scenes > 1){
      #     print("More than one scene in image file.")
      #   }
      # }
      #
      # # Check for phases (with apotome)
      # if("H" %in% axes){
      #   pos_phases <- grep(pattern = "H", x = axes)
      #   number_phases <- dim(image_loaded)[pos_phases]
      #
      #   if(!apotome){
      #     print("Apotome was used. Parameter is set TRUE.")
      #     apotome <- TRUE
      #   }
      # }

      # Work with apotome image ###
      # Reduce apotome layers to one
      # Algorithm taken from SCHAEFER et al. (2004)
      # ["Structured illumination microscopy: artefact analysis and
      #   reduction utilizing a parameter optimization approach"]

      # if(apotome){
      #
      #   if(apotome_section){
      #     image_sectioned <- image_loaded[1,,,,,]
      #   }else{
      #     image_conventional <- image_loaded[1,,,,,]
      #   }
      #
      #
      #   for(chan in 1:number_of_channels){
      #
      #     cos_term <- 0
      #     sin_term <- 0
      #     conv_term <- 0
      #
      #     for(phase in 1:number_phases){
      #
      #       if(pos_phases == 1 & pos_channels == 3){
      #         if(apotome_section){
      #           # cos terms
      #           cos_term <- cos_term + image_loaded[phase,,chan,,,] *
      #             cos(2*pi*(phase-1)/number_phases)
      #
      #           # sin terms
      #           sin_term <- sin_term + image_loaded[phase,,chan,,,] *
      #             sin(2*pi*(phase-1)/number_phases)
      #         }else{
      #           # conventional term
      #           conv_term <- conv_term + image_loaded[phase,,chan,,,]
      #         }
      #
      #       }else{
      #         print("Position of phases not 1.")
      #       }
      #
      #     }
      #
      #     cos_term <- (2/number_phases) * cos_term
      #     sin_term <- (2/number_phases) * sin_term
      #     conv_term <- conv_term / number_phases
      #
      #     if(apotome_section){
      #       image_sectioned[chan,,] <- sqrt(cos_term^2+sin_term^2)
      #     }else{
      #       image_conventional[chan,,] <- conv_term
      #     }
      #
      #   }
      #
      #   if(apotome_section){
      #     image_loaded <- image_sectioned
      #   }else{
      #     image_loaded <- image_conventional
      #   }
      #
      # }

      # Get bit depth of camera
      # camera_bit_depth <- gsub(
      #   pattern = ".+<ComponentBitCount>(.+)</ComponentBitCount>.+",
      #   replacement = "\\1",
      #   x = metadata)
      # camera_bit_depth <- as.numeric(camera_bit_depth)
      # if(!is.na(camera_bit_depth) && (camera_bit_depth != 0)){
      #   bit_depth <- camera_bit_depth
      # }

      # Find red, green, and blue channel ID
      # wavelengths <- gsub(
      #   pattern =  paste(".+<EmissionWavelength>(.+)</EmissionWavelength>.+",
      #                    ".+<EmissionWavelength>(.+)</EmissionWavelength>.+",
      #                    ".+<EmissionWavelength>(.+)</EmissionWavelength>.+",
      #                    sep=""),
      #   replacement = "\\1,\\2,\\3",
      #   x = metadata)
      # wavelengths <- as.numeric(strsplit(wavelengths, split = ",")[[1]])
      #
      # red_id <- which(wavelengths == max(wavelengths))
      # blue_id <- which(wavelengths == min(wavelengths))
      #
      # if(length(wavelengths) == 3){
      #   green_id <- c(1:3)[!(c(1:3) %in% red_id | c(1:3) %in% blue_id)]
      # }else{
      #   print("There are more or fewer than 3 wavelengths.")
      #   return()
      # }
      #
      # rgb_layers <- c(red_id, green_id, blue_id)
      #
      # rm(czi_class)

      # Make a 3-D array of the image

      # if(number_of_channels == 3){
      #
      #   # Delete the dimensions of an array which have only one level
      #   image_loaded <- drop(image_loaded)
      #
      #
      #   # New position of X,Y,Channels
      #   pos_channels <- which(dim(image_loaded) == number_of_channels)
      #   pos_x <- which(dim(image_loaded) == dim_x)
      #   pos_y <- which(dim(image_loaded) == dim_y)
      #
      #   ## Permute dimensions of array
      #   #pos_x <- pos_x - min(pos_x, pos_y, pos_channels) + 1
      #   #pos_y <- pos_y - min(pos_x, pos_y, pos_channels) + 1
      #   #pos_channels <- pos_channels - min(pos_x, pos_y, pos_channels) + 1
      #
      #   image_loaded <- aperm(a = image_loaded, c(pos_y, pos_x, pos_channels))
      #
      #   # Reorder the layers according to the colors
      #   copy_image_loaded <- image_loaded
      #
      #   image_loaded[,,1] <- copy_image_loaded[,,rgb_layers[1]]
      #   image_loaded[,,2] <- copy_image_loaded[,,rgb_layers[2]]
      #   image_loaded[,,3] <- copy_image_loaded[,,rgb_layers[3]]
      #
      #   rm(copy_image_loaded)
      #
      # }else{
      #   print(paste("We do not have a tif-image with three layers",
      #               " representing red, green, and blue.", sep=""))
      #   return()
      # }

      # image_loaded <- image_loaded/(2^bit_depth-1)


      # Get the exposure time for every layer
      # exposure_time <- gsub(
      #   pattern =  paste(".+<ExposureTime>(.+)</ExposureTime>.+",
      #                    ".+<ExposureTime>(.+)</ExposureTime>.+",
      #                    ".+<ExposureTime>(.+)</ExposureTime>.+",
      #                    sep=""),
      #   replacement = "\\1,\\2,\\3",
      #   x = metadata)
      #
      # exposure_time <- unlist(strsplit(x = exposure_time, split = ","))

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



      # Reorder the channels depending on the channel information


      # if(is.list(image_loaded)){
      #   # Combine channels into one array
      #   image_dummy <- array(dim = c(dim(image_loaded[[1]])[1],
      #                                dim(image_loaded[[1]])[2],
      #                                length(image_loaded) ), data = 0)
      #
      #   for(layer in 1:length(image_loaded)){
      #     image_dummy[,,layer] <- image_loaded[[layer]]
      #   }
      #
      #   image_loaded <- image_dummy
      #   rm(image_dummy)
      #
      # }

      # Load metadata


    }

    # Reduce image size if dim_x % 4 != 0 or dim_y % 4 != 0
    # (The implementation of EBImage::clahe assumes that the X- and Y image
    # dimensions are an integer multiple of the X- and Y sizes of the
    # contextual regions)

    number_of_contextual_regions <- 4

    # dim x
    number_of_pixels_to_crop_x <- dim(image_loaded)[1] %% number_of_contextual_regions
    number_of_pixels_to_crop_y <- dim(image_loaded)[2] %% number_of_contextual_regions
    image_cropped <- "no"

    if(number_of_pixels_to_crop_x != 0 || number_of_pixels_to_crop_y != 0){
      print(paste0("The x or y dimension is not a multiple of ",
                   number_of_contextual_regions,
                   ". Therefore, we are cropping the image."))

      image_cropped <- "yes"

      if(length(dim(image_loaded)==3)){
        image_loaded <- image_loaded[
          1:(dim(image_loaded)[1]-number_of_pixels_to_crop_x),
          1:(dim(image_loaded)[2]-number_of_pixels_to_crop_y),]
      }else if(length(dim(image_loaded)==2)){
        image_loaded <- image_loaded[
          1:(dim(image_loaded)[1]-number_of_pixels_to_crop_x),
          1:(dim(image_loaded)[2]-number_of_pixels_to_crop_y)]
      }

    }


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
    #     print(paste("Bit depth is higher than 16. Something might",
    #                 " be wrong.", sep=""))
    #     return()
    #   }
    #
    # }


    # Number of total pixels
    number_of_total_pixels <- dim(image_loaded)[1]*dim(image_loaded)[2]

    # -------------------------------------------------------------------- #
    # ------------------ Find background intensity ----------------------- #
    # -------------------------------------------------------------------- #

    # Go through every layer and save the foreground as a mask
    number_of_channels <- df_metadata$number_of_channels[1]

    for(j in 1:number_of_channels){
      image_dummy <- image_loaded[,,j]
      image_dummy <- EBImage::gblur(image_dummy, sigma = 5)

      if(j == 1){
        mask_foreground <- EBImage::thresh(image_dummy, w=50, h=50, offset=0.001)
      }else{

        # Add the foreground to one mask
        mask_foreground_dummy <- EBImage::thresh(image_dummy, w=50, h=50, offset=0.001)
        mask_foreground <- mask_foreground + as.numeric(mask_foreground_dummy > mask_foreground)
      }

    }
    rm(j)
    rm(image_dummy)
    rm(mask_foreground_dummy)

    # Number of pixels of the foreground
    number_of_pixels_foreground <- sum(mask_foreground)

    # -------------------------------------------------------------------- #
    # ---------------------- Data manipulation --------------------------- #
    # -------------------------------------------------------------------- #

    # Normalize intensity --------------------------------------------------
    image_normalized <- readCzi::normalizeIntensity(image = image_loaded)

    # Get Background of normalized image and loaded image
    image_normalized_foreground <- image_normalized
    image_normalized_background <- image_normalized
    image_background <- image_loaded
    image_foreground <- image_loaded

    for(j in 1:number_of_channels){
      image_normalized_foreground[,,j] <- image_normalized[,,j]*mask_foreground
      image_normalized_background[,,j] <- image_normalized[,,j]*(-1*mask_foreground+1)
      image_background[,,j] <- image_loaded[,,j]*(-1*mask_foreground+1)
      image_foreground[,,j] <- image_loaded[,,j]*(mask_foreground)
    }
    rm(j)

    # Use Contrast Limited Adaptive Histogram Equalization
    image_histogram_equalization <- EBImage::clahe(x = image_loaded, nx = number_of_contextual_regions)

    # Find the nuclei ------------------------------------------------------

    # Save only color layer of nuclei
    if(use_histogram_equalized){
      image_nuclei <- getLayer(image = image_histogram_equalization, layer = nucleus_color)
    }else{
      image_nuclei <- getLayer(image = image_loaded, layer = nucleus_color)
    }

    # Brighten nuclei image for apotome section image
    if(normalize_nuclei_layer || apotome_section){
      image_nuclei <- image_nuclei/max(image_nuclei)
    }

    Image_nuclei <- EBImage::Image(image_nuclei)
    rm(image_nuclei)
    #display(Image_nuclei)

    # Blur the image
    if(is.null(blur_sigma)){
      if(apotome_section){
        Image_nuclei <- EBImage::gblur(Image_nuclei, sigma = 14)
      }else{
        # Image_nuclei <- EBImage::gblur(Image_nuclei, sigma = 7)
        Image_nuclei <- EBImage::gblur(Image_nuclei, sigma = 6)
      }
    }else if(blur_sigma > 0){
      Image_nuclei <- EBImage::gblur(Image_nuclei, sigma = blur_sigma)
    }


    if(is.null(magnification_objective)){
      magnification <- df_metadata$objective_magnification[df_metadata$fileName == file_names[i]]
    }else{
      magnification <- magnification_objective
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
        nmask <- EBImage::thresh(Image_nuclei, w=500, h=500, offset=0.03)
      }
    }else{
      nmask <- EBImage::thresh(Image_nuclei, w=thresh_w_h_nuc, h=thresh_w_h_nuc, offset=thresh_offset)
    }
    #display(nmask)


    # Morphological opening to remove objects smaller than the structuring element
    nmask <- EBImage::opening(nmask, EBImage::makeBrush(5, shape='disc'))
    # Fill holes
    nmask <- EBImage::fillHull(nmask)

    # Save black-and-white-image as nucleus mask
    nucleus_mask <- nmask

    # Number of pixels of the nucleus mask
    number_of_pixels_nucleus_region <- sum(nucleus_mask)

    # Number of pixels of the foreground without the nucleus part
    mask_foreground_wo_nucleus <- mask_foreground - nucleus_mask
    mask_foreground_wo_nucleus[mask_foreground_wo_nucleus < 0] <- 0

    number_of_pixels_foreground_without_nucleus_region <- sum(mask_foreground_wo_nucleus)

    # Label each connected set of pixels with a distinct ID
    nmask <- EBImage::bwlabel(nmask)

    #display(nmask)


    # Watershed in order to distinct nuclei that are too close to each other
    # nmask_watershed <-  EBImage::watershed(
    #   EBImage::distmap(nmask), tolerance = 0.4, ext = 2)

    # if(magnification < 20){
    #   nmask_watershed <-  EBImage::watershed(
    #     EBImage::distmap(nmask), tolerance = 0.1, ext = 3)
    # }else if(magnification < 40){
    #   nmask_watershed <-  EBImage::watershed(
    #     EBImage::distmap(nmask), tolerance = 0.1, ext = 6)
    # }else{
    #   nmask_watershed <-  EBImage::watershed(
    #     EBImage::distmap(nmask), tolerance = 0.1, ext = 9)
    # }
    #
    # if(magnification < 20){
    #   nmask_watershed <-  EBImage::watershed(
    #     EBImage::distmap(nmask), tolerance = 0.45, ext = 3) # it was tolerance = 0.5
    # }else if(magnification < 40){
    #   nmask_watershed <-  EBImage::watershed(
    #     EBImage::distmap(nmask), tolerance = 0.45, ext = 2)
    # }else{
    #   nmask_watershed <-  EBImage::watershed(
    #     EBImage::distmap(nmask), tolerance = 0.3, ext = 3) # used to be 045, 9
    # }

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
        EBImage::distmap(nmask), tolerance = 0.1, ext = 6) # used to be 045, 9 | 0.3, 3

      # display(colorLabels(EBImage::watershed(EBImage::distmap(nmask), tolerance = 0.1, ext = 6)), all=TRUE)
    }


    #display(colorLabels(nmask_watershed), all=TRUE)
    #display(colorLabels(nmask), all=TRUE)
    # display(colorLabels(EBImage::watershed(EBImage::distmap(nmask), tolerance = 0.3, ext = 6)), all=TRUE)

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
    nucNo <- max(nmask_watershed)

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
          number_color = "green")
        Image_nuclei_numbers <- addNumberToImage(
          image = Image_nuclei_numbers,
          number = df_nmask_watershed$newNucNo[j],
          pos_x = df_nmask_watershed$x_coord[j],
          pos_y = df_nmask_watershed$y_coord[j],
          number_size_factor = number_size_factor,
          number_color = "blue")
      }
      rm(j)

    }

    mean_nuc_size <- mean(df_nmask_watershed$Freq)


    # if(length(table_nmask_watershed[-1]) > 0){
    #   # remove 0
    #   nuc_numbers <- as.integer(names(table_nmask_watershed[-1]))

    # for(j in 1:length(nuc_numbers)){
    #
    #   # Find approximate midpoint of every nucleus
    #   dummy_coordinates <- which(
    #     EBImage::imageData(nmask_watershed) ==
    #       nuc_numbers[j], arr.ind = TRUE)
    #
    #
    #   pos_x <- round(mean(dummy_coordinates[,1]))
    #   pos_y <- round(mean(dummy_coordinates[,2]))
    #
    #   Image_nuclei_numbers <- addNumberToImage(
    #     image = Image_nuclei_numbers,
    #     number = j,
    #     pos_x = pos_x,
    #     pos_y = pos_y,
    #     number_size_factor = number_size_factor,
    #     number_color = "green")
    #   Image_nuclei_numbers <- addNumberToImage(
    #     image = Image_nuclei_numbers,
    #     number = j,
    #     pos_x = pos_x,
    #     pos_y = pos_y,
    #     number_size_factor = number_size_factor,
    #     number_color = "blue")
    # }
    # rm(j)
    # }

    # Add border of nuclei and save file
    Image_nuclei <- EBImage::Image(Image_nuclei, colormode = "color")
    # EBImage::colorMode(Image_nuclei) <- "color"
    # display(Image_nuclei)

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

      if(use_histogram_equalized){
        image_protein_in_nuc <- EBImage::clahe(x = image_protein_in_nuc, nx = 4)
        #display(image_protein_in_nuc)
      }

      Image_protein_in_nuc <- EBImage::Image(image_protein_in_nuc)
      rm(image_protein_in_nuc)
      #display(Image_protein_in_nuc)

      # Blur the image
      if(is.null(blur_sigma)){
        if(apotome_section){
          Image_protein_in_nuc <- EBImage::gblur(Image_protein_in_nuc, sigma = 14)
        }else{
          Image_protein_in_nuc <- EBImage::gblur(Image_protein_in_nuc, sigma = 6)
        }
      }else if(blur_sigma > 0){
        Image_protein_in_nuc <- EBImage::gblur(Image_protein_in_nuc, sigma = blur_sigma)
      }

      #display(Image_protein_in_nuc)

      # Mask the proteins within the nucleus
      # Todo: aufteilen nach beiden Bedingungen und dann Fehlendes setzen
      if(is.null(thresh_w_h_nuc)){
        if(magnification < 20){
          thresh_w_h_nuc <- 15
        }else if(magnification < 40){
          thresh_w_h_nuc <- 30
        }else if(magnification < 60){
          thresh_w_h_nuc <- 60
        }else{
          thresh_w_h_nuc <- 500
        }
      }

      if(is.null(thresh_offset_protein_in_nucleus)){
        thresh_offset_protein_in_nucleus <- mean(Image_protein_in_nuc)+0.1
      }

      pmask <- EBImage::thresh(Image_protein_in_nuc,
                               w=thresh_w_h_nuc,
                               h=thresh_w_h_nuc,
                               offset=thresh_offset_protein_in_nucleus)

      #display(pmask)

      # Morphological opening to remove objects smaller than the structuring element
      pmask <- EBImage::opening(pmask, EBImage::makeBrush(5, shape='disc'))
      # Fill holes
      pmask <- EBImage::fillHull(pmask)
      #display(pmask)


      # Combine nmask_watershed and pmask and count the resulting nuclei containing
      # the proteins we are looking for
      n_p_mask <- nmask_watershed*pmask
      #display(n_p_mask)

      # remove objects that are smaller than 0.1*nuc_size
      table_npmask <- table(n_p_mask)
      # table_nmask_watershed
      to_be_removed <- as.integer(names(
        which(table_npmask < 0.01 * table_nmask[match(names(table_npmask), names(table_nmask))])))
      if(length(to_be_removed) > 0){
        n_p_mask[n_p_mask %in% to_be_removed] <- 0
      }

      # display(n_p_mask)


      # Combine possibly multiple positive nuclei
      if(magnification < 20){
        n_p_mask_watershed <-  EBImage::watershed(
          EBImage::distmap(n_p_mask), tolerance = 0.45, ext = 3) # it was tolerance = 0.5
      }else if(magnification < 40){
        n_p_mask_watershed <-  EBImage::watershed(
          EBImage::distmap(n_p_mask), tolerance = 0.1, ext = 6) #sqrt(15), see thresh(Image_nuclei...
      }else if(magnification < 60){
        n_p_mask_watershed <-  EBImage::watershed(
          EBImage::distmap(n_p_mask), tolerance = 0.1, ext = 6) # used to be 045, 9 | 0.3, 3
      }else{
        n_p_mask_watershed <-  EBImage::watershed(
          EBImage::distmap(n_p_mask), tolerance = 0.1, ext = 6) # used to be 045, 9 | 0.3, 3
      }

      # display(colorLabels(n_p_mask_watershed))


      # remove objects that are smaller than 0.1*min_nuc_size
      # table_npmask <- table(n_p_mask)
      table_npmask <- table(n_p_mask_watershed)





      # Count the nuclei containing the proteins
      # positive_nuclei <- as.numeric(names(table(n_p_mask)))
      positive_nuclei <- as.numeric(names(table(n_p_mask_watershed)))
      positive_nuclei <- positive_nuclei[positive_nuclei != 0]
      # positive_nuclei <- sort(df_nmask_watershed$newNucNo[match(positive_nuclei, df_nmask_watershed$NucNo)])

      # Attention: We have a new list now. The numbers do not correspondend
      # to the actual nuclei numbers
      # (We could no map both the nuclei and the positive nuclei map and
      # only keep the number of the nuclei with the largest coverage.)

      # print(positive_nuclei)
      # n_p_mask <- EBImage::bwlabel(n_p_mask)
      #display(n_p_mask)
      # Count number of cells
      # nuc_with_proteins_No <- max(n_p_mask)

      nuc_with_proteins_No <- length(positive_nuclei)

      # Add border of nuclei with proteins and save file
      Image_nuclei_numbers_proteins <- Image_nuclei_numbers
      EBImage::colorMode(Image_nuclei_numbers_proteins) <- "color"

      Image_nuclei_numbers_proteins <- EBImage::paintObjects(
        # x = n_p_mask,
        x = n_p_mask_watershed,
        tgt = Image_nuclei_numbers_proteins,
        opac = c(1,0.8),
        col=c('#D7FE14','#D7FE14'))
      # display(Image_nuclei_numbers_proteins)

      # Display the number of nuclei with proteins
      print(paste("Number of nuclei that contain other colored proteins: ",
                  nuc_with_proteins_No, sep=""))

    }

    # Count the number of cells that contain a third colored protein  ------

    if(!is.null(protein_in_cytosol_color) & protein_in_cytosol_color != "none"){


      # Save only color layer of the second protein colored
      image_protein_in_cytosol <- getLayer(
        image = image_loaded, layer = protein_in_cytosol_color)

      # Brighten cytosol image for apotome section image
      if(apotome_section){
        image_protein_in_cytosol <- image_protein_in_cytosol/max(image_protein_in_cytosol)
      }

      Image_protein_in_cytosol <- EBImage::Image(image_protein_in_cytosol)
      rm(image_protein_in_cytosol)
      #display(Image_protein_in_cytosol)

      # Blur the image
      if(apotome_section){
        Image_protein_in_cytosol <- EBImage::gblur(Image_protein_in_cytosol, sigma = 10)
      }else{
        Image_protein_in_cytosol <- EBImage::gblur(Image_protein_in_cytosol, sigma = 5)
      }

      #display(Image_protein_in_cytosol)

      # Mask the proteins within the cytosol
      # if(grepl(pattern = "_20x_", file_names[i])){
      #   # Smaller moving rectangle if the magnification is 20x (instead of 40x)
      #   cytosolmask <- EBImage::thresh(Image_protein_in_cytosol, w=100, h=100, offset=0.01)
      # }else{
      #   cytosolmask <- EBImage::thresh(Image_protein_in_cytosol, w=200, h=200, offset=0.01)
      # }

      mean_intensity <- mean(Image_protein_in_cytosol)

      if(is.null(thresh_w_h_nuc) || is.null(thresh_offset_protein_in_cytosol)){
        if(magnification < 20){
          cytosolmask <- EBImage::thresh(Image_protein_in_cytosol, w=50, h=50, offset=2*mean_intensity) #0.2
        }else if(magnification < 40){
          # Smaller moving rectangle if the magnification is 20x (instead of 40x)
          cytosolmask <- EBImage::thresh(Image_protein_in_cytosol, w=100, h=100, offset=2*mean_intensity)
        }else{
          cytosolmask <- EBImage::thresh(Image_protein_in_cytosol, w=200, h=200, offset=2*mean_intensity)
        }
      }else{
        cytosolmask <- EBImage::thresh(Image_protein_in_cytosol,
                                       w=thresh_w_h_nuc,
                                       h=thresh_w_h_nuc,
                                       offset=thresh_offset_protein_in_cytosol)
      }

      #display(cytosolmask)

      # Morphological opening to remove objects smaller than the structuring element
      cytosolmask <- EBImage::opening(cytosolmask, EBImage::makeBrush(5, shape='disc'))
      # Fill holes
      cytosolmask <- EBImage::fillHull(cytosolmask)
      #display(cytosolmask)

      # Keep only those cell bodies that contain a nucleus
      # cytosolmask <- EBImage::bwlabel(cytosolmask)
      # no_of_cytosols <- max(cytosolmask)

      # Go through every stained cytosol and check for nucleus
      n_c_mask <- cytosolmask * nmask_watershed

      # remove objects that are smaller than 0.1*min_nuc_size
      table_ncmask <- table(n_c_mask)
      to_be_removed <- as.integer(names(which(table_ncmask < 0.1*nuc_min_size)))

      n_c_mask[n_c_mask %in% to_be_removed] <- 0
      # display(n_c_mask)

      # Count the cells containing the proteins
      positive_cells <- as.numeric(names(table(n_c_mask)))
      positive_cells <- positive_cells[positive_cells != 0]
      positive_cells <- sort(df_nmask_watershed$newNucNo[match(positive_cells, df_nmask_watershed$NucNo)])
      # print(positive_cells)
      # n_p_mask <- EBImage::bwlabel(n_p_mask)

      cell_with_proteins_No <- length(positive_cells)

      # for(j in 1:no_of_cytosols){
      #   collocation_found <- max((cytosolmask==j)*nmask)
      #   if(collocation_found == 0){
      #     cytosolmask[cytosolmask==j] <- 0
      #   }
      #   rm(j)
      # }

      # cytosolmask <- EBImage::bwlabel(cytosolmask)
      # no_of_cytosols <- max(cytosolmask)

      # Combine nmask and cytosolmask and count the resulting cell bodies
      # (nuclei that are within the stained proteins in the cytosol)
      # n_c_mask <- nmask_watershed*(!cytosolmask==0)
      #display(n_c_mask)

      # table_n_c_mask <- table(n_c_mask)

      # Count number of cells containing staining for protein
      # cell_with_proteins_No <- length(names(table_n_c_mask))-1
      #display(n_c_mask)


      # remove objects that are smaller than min_nuc_size
      # to_be_removed <- as.integer(names(which(table_n_c_mask < nuc_min_size)))
      # if(length(to_be_removed) > 0){
      #   for(j in 1:length(to_be_removed)){
      #     EBImage::imageData(n_c_mask)[
      #       EBImage::imageData(n_c_mask) == to_be_removed[j]] <- 0
      #   }
      #   rm(j)
      # }

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

    # Detect clusters of membrane ligands  ---------------------------------

    if(protein_in_membrane_color != "none" & apotome_section){

      # Save only color layer membrane ligands
      if(use_histogram_equalized){
        image_membrane_ligands <- getLayer(image = image_histogram_equalization, layer = protein_in_membrane_color)
      }else{
        image_membrane_ligands <- getLayer(image = image_loaded, layer = protein_in_membrane_color)
      }

      # Normalize membrane ligands image
      image_membrane_ligands <- image_membrane_ligands/max(image_membrane_ligands)

      Image_membrane_ligands <- EBImage::Image(image_membrane_ligands)
      rm(image_membrane_ligands)
      #display(Image_membrane_ligands)



      # # Blur the image
      # if(is.null(blur_sigma)){
      #   if(apotome_section){
      #     Image_nuclei <- EBImage::gblur(Image_nuclei, sigma = 14)
      #   }else{
      #     Image_nuclei <- EBImage::gblur(Image_nuclei, sigma = 7)
      #   }
      # }else if(blur_sigma > 0){
      #   Image_nuclei <- EBImage::gblur(Image_nuclei, sigma = blur_sigma)
      # }


      # Mask the clusters
      # if(grepl(pattern = "_63x_", file_names[i])){
      if(magnification == 63){
        cmask <- EBImage::thresh(Image_membrane_ligands, w=500, h=500, offset=0.19)
      }else{
        print("magnification of objective other than 63x used")
        cmask <- EBImage::thresh(Image_membrane_ligands, w=10, h=10, offset=0.3)
      }

      #display(cmask)

      # Fill holes
      cmask <- EBImage::fillHull(cmask)

      # Save black-and-white-image as cluster mask
      cluster_mask <- cmask
      #
      #       # Number of pixels of the nucleus mask
      #       number_of_pixels_nucleus_region <- sum(nucleus_mask)
      #
      #       # Number of pixels of the foreground without the nucleus part
      #       mask_foreground_wo_nucleus <- mask_foreground - nucleus_mask
      #       mask_foreground_wo_nucleus[mask_foreground_wo_nucleus < 0] <- 0
      #
      #       number_of_pixels_foreground_without_nucleus_region <- sum(mask_foreground_wo_nucleus)

      # Label each connected set of pixels with a distinct ID
      cmask <- EBImage::bwlabel(cmask)

      #display(cmask)

      # # Record all nuclei that are at the edges of the image
      # left  <- table(nmask[1:number_of_pixels_at_border_to_disregard,
      #                      1:dim(nmask)[2]])
      # top   <- table(nmask[1:dim(nmask)[1],
      #                      1:number_of_pixels_at_border_to_disregard])
      # right <- table(nmask[
      #   (dim(nmask)[1]-number_of_pixels_at_border_to_disregard):dim(nmask)[1],
      #   1:dim(nmask)[2]])
      # bottom <- table(nmask[
      #   1:dim(nmask)[1],
      #   (dim(nmask)[2]-number_of_pixels_at_border_to_disregard):dim(nmask)[2]])
      #
      #
      # left <- as.integer(names(left))
      # top <- as.integer(names(top))
      # right <- as.integer(names(right))
      # bottom <- as.integer(names(bottom))
      #
      # nuclei_at_borders <- unique(c(left, top, right, bottom))
      #
      # # delete the 0 in the list (these are pixels that do not belong to a nucleus)
      # nuclei_at_borders <- nuclei_at_borders[nuclei_at_borders != 0]
      #
      # # Delete all nuclei at border
      # if(length(nuclei_at_borders) > 0){
      #   for(j in 1:length(nuclei_at_borders)){
      #     EBImage::imageData(nmask)[
      #       EBImage::imageData(nmask) == nuclei_at_borders[j]] <- 0
      #   }
      #   rm(j)
      # }

      #display(nmask)

      # Delete all remaining nuclei that are smaller than 5% of the median size
      # object sizes
      # barplot(table(nmask)[-1])

      # table_nmask <- table(nmask)
      # nuc_min_size <- 0.1*stats::median(table_nmask[-1])
      #
      # # remove objects that are smaller than min_nuc_size
      # to_be_removed <- as.integer(names(which(table_nmask < nuc_min_size)))
      # if(length(to_be_removed) > 0){
      #   for(j in 1:length(to_be_removed)){
      #     EBImage::imageData(nmask)[
      #       EBImage::imageData(nmask) == to_be_removed[j]] <- 0
      #   }
      #   rm(j)
      # }
      #
      # # Recount nuclei
      # nmask <- EBImage::bwlabel(nmask)
      # #display(nmask)

      # Watershed in order to distinct clusters that are too close to each other
      cmask_watershed <-  EBImage::watershed(
        EBImage::distmap(cmask), tolerance = 0.1, ext = 1)
      #display(colorLabels(cmask_watershed), all=TRUE)

      # Count number of cluster
      clustersNo <- max(cmask_watershed)

      # Mean size of clusters
      cluster_sizes <- as.vector(table(cmask)[-1])
      mean_cluster_size <- mean(cluster_sizes)
      median_cluster_size <- median(cluster_sizes)

      # TODO: save image with and without the cluster mask

      # Use the cluster mask to leave only the clusters part of the images
      Image_clusters_part <- EBImage::Image(image_normalized) * EBImage::toRGB(cluster_mask)

      # Use the nucleus mask to cut out nuclei of the images
      cluster_mask <- 1-cluster_mask
      Image_non_clusters_part <- EBImage::Image(image_normalized) * EBImage::toRGB(cluster_mask)

    }else{
      clustersNo <- 0
      mean_cluster_size <- NA
      median_cluster_size <- NA
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
      df_results[i,"color_of_third_protein_in_cytosol"] <- protein_in_cytosol_color
      df_results[i,"number_of_cells_with_third_protein"] <- cell_with_proteins_No
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

    df_results[i,"intensity_sum_red_foreground"] <- sum((Image_loaded * EBImage::toRGB(mask_foreground))[,,1])
    df_results[i,"intensity_sum_green_foreground"] <- sum((Image_loaded * EBImage::toRGB(mask_foreground))[,,2])
    df_results[i,"intensity_sum_blue_foreground"] <- sum((Image_loaded * EBImage::toRGB(mask_foreground))[,,3])

    df_results[i,"intensity_sum_red_foreground_without_nucleus_region"] <- sum((Image_loaded * EBImage::toRGB(mask_foreground_wo_nucleus))[,,1])
    df_results[i,"intensity_sum_green_foreground_without_nucleus_region"] <- sum((Image_loaded * EBImage::toRGB(mask_foreground_wo_nucleus))[,,2])
    df_results[i,"intensity_sum_blue_foreground_without_nucleus_region"] <- sum((Image_loaded * EBImage::toRGB(mask_foreground_wo_nucleus))[,,3])

    df_results[i,"intensity_mean_red_background"] <- mean(image_background[,,1])
    df_results[i,"intensity_mean_green_background"] <- mean(image_background[,,2])
    df_results[i,"intensity_mean_blue_background"] <- mean(image_background[,,3])

    df_results[i,"number_of_total_pixels"] <- number_of_total_pixels
    df_results[i,"number_of_pixels_nucleus_region"] <- number_of_pixels_nucleus_region
    df_results[i,"number_of_pixels_foreground"] <-  number_of_pixels_foreground
    df_results[i,"number_of_pixels_foreground_without_nucleus_region"] <- number_of_pixels_foreground_without_nucleus_region

    df_results[i, "number_of_clusters"] <- clustersNo
    df_results[i, "mean_cluster_size"] <- mean_cluster_size
    df_results[i, "median_cluster_size"] <- median_cluster_size

    df_results[i, "image_cropped"] <- image_cropped


    # if(image_format == "czi"){
    #   if("laser_scan_pixel_time_1" %in% names(df_metadata)){
    #     df_results[i,"exposure_or_laser_scan_pixel_time_channel_1"] <- df_metadata$laser_scan_pixel_time_1[1]
    #     df_results[i,"exposure_or_laser_scan_pixel_time_channel_2"] <- df_metadata$laser_scan_pixel_time_1[2]
    #     df_results[i,"exposure_or_laser_scan_pixel_time_channel_3"] <- df_metadata$laser_scan_pixel_time_1[3]
    #   }else if("exposure_time_1" %in% names(df_metadata)){
    #     df_results[i,"exposure_or_laser_scan_pixel_time_channel_1"] <- df_metadata$exposure_time_1[1]
    #     df_results[i,"exposure_or_laser_scan_pixel_time_channel_2"] <- df_metadata$exposure_time_2[2]
    #     df_results[i,"exposure_or_laser_scan_pixel_time_channel_3"] <- df_metadata$exposure_time_3[3]
    #   }
    #
    # }

    # Save all images ------------------------------------------------------

    # if(image_format == "czi"){
      # Get information for adding a scale bar
      # length_per_pixel_x <- gsub(
      #   pattern = paste(".+<Items>[[:space:]]+<Distance Id=\"X\">[[:space:]]+",
      #                   "<Value>(.{2,25})</Value>.+",sep=""),
      #   replacement = "\\1",
      #   x = metadata)
      # length_per_pixel_x <- tolower(length_per_pixel_x)
      # length_per_pixel_x <- as.numeric(length_per_pixel_x)
      length_per_pixel_x_in_um <- df_metadata$scaling_x_in_um[1]

      # length_per_pixel_y <- gsub(
      #   pattern = paste(".+<Items>.+<Distance Id=\"Y\">[[:space:]]+",
      #                   "<Value>(.{2,25})</Value>.+",sep=""),
      #   replacement = "\\1",
      #   x = metadata)
      # length_per_pixel_y <- tolower(length_per_pixel_y)
      # length_per_pixel_y <- as.numeric(length_per_pixel_y)
      length_per_pixel_y_in_um <- df_metadata$scaling_y_in_um[1]

      if(length_per_pixel_x_in_um != length_per_pixel_y_in_um){
        print("Dimension in x- and y-directions are different! ERROR!")
      }

      # # Save metadata in txt file
      # utils::write.table(metadata, file = paste(output_dir,image_name_wo_czi,
      #                                           "_metadata.txt", sep = ""),
      #                    sep = "", row.names = FALSE, col.names = FALSE, quote = FALSE)


      # Original image (converted to tif)
      if(add_scale_bar){
        image_loaded <- addScaleBar(image = image_loaded,
                                    length_per_pixel = length_per_pixel_x_in_um)
      }
    # }

    # tiff::writeTIFF(what = image_loaded,
    #                 where = paste(output_dir,
    #                               image_name_wo_czi,
    #                               "_original.tif",
    #                               sep = ""),
    #                 bits.per.sample = 8L, compression = "none",
    #                 reduce = TRUE)

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

    # tiff::writeTIFF(what = image_normalized,
    #                 where = paste(output_dir,
    #                               image_name_wo_czi,
    #                               "_normalized.tif",
    #                               sep = ""),
    #                 bits.per.sample = 8L, compression = "none",
    #                 reduce = TRUE)
    Image_normalized <- EBImage::Image(data = image_normalized, colormode = "Color")
    EBImage::writeImage(x = Image_normalized,
                        files = paste(output_dir, image_name_wo_czi, "_normalized.tif", sep = ""),
                        type = "tiff",
                        bits.per.sample = 8)

    # tiff::writeTIFF(what = image_normalized_background,
    #                 where = paste(output_dir,
    #                               image_name_wo_czi,
    #                               "_normalized_background.tif",
    #                               sep = ""),
    #                 bits.per.sample = 8L, compression = "none",
    #                 reduce = TRUE)
    Image_normalized_background <- EBImage::Image(data = image_normalized_background, colormode = "Color")
    EBImage::writeImage(x = Image_normalized_background,
                        files = paste(output_dir, image_name_wo_czi, "_normalized_background.tif", sep = ""),
                        type = "tiff",
                        bits.per.sample = 8)

    # tiff::writeTIFF(what = image_histogram_equalization,
    #                 where = paste(output_dir,
    #                               image_name_wo_czi,
    #                               "_histogram_equalized.tif",
    #                               sep = ""),
    #                 bits.per.sample = 8L, compression = "none",
    #                 reduce = TRUE)
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
    # tiff::writeTIFF(what = Image_nuclei,
    #                 where = paste(output_dir,
    #                               image_name_wo_czi,
    #                               "_nuclei.tif",
    #                               sep = ""),
    #                 bits.per.sample = 8L, compression = "none",
    #                 reduce = TRUE)
    # Image_nuclei <- EBImage::Image(data = Image_nuclei, colormode = "Color")
    EBImage::writeImage(x = Image_nuclei,
                        files = paste(output_dir, image_name_wo_czi, "_nuclei.tif", sep = ""),
                        type = "tiff",
                        bits.per.sample = 8)

    # tiff::writeTIFF(what = Image_nuclei_numbers,
    #                 where = paste(output_dir,
    #                               image_name_wo_czi,
    #                               "_nuclei_numbers.tif",
    #                               sep = ""),
    #                 bits.per.sample = 8L, compression = "none",
    #                 reduce = TRUE)
    # Image_nuclei_numbers <- EBImage::Image(data = Image_nuclei_numbers, colormode = "Color")
    EBImage::writeImage(x = Image_nuclei_numbers,
                        files = paste(output_dir, image_name_wo_czi, "_nuclei_numbers.tif", sep = ""),
                        type = "tiff",
                        bits.per.sample = 8)

    # Images with marked nuclei and borders around the second protein in nuc
    if(add_scale_bar){
      Image_nuclei_numbers_proteins <- addScaleBar(
        image = Image_nuclei_numbers_proteins,
        length_per_pixel = length_per_pixel_x_in_um)
    }
    if(!is.null(protein_in_nuc_color) & protein_in_nuc_color != "none"){
      # tiff::writeTIFF(what = Image_nuclei_numbers_proteins,
      #                 where = paste(output_dir,
      #                               image_name_wo_czi,
      #                               "_nuclei_numbers_proteins1_nuc.tif",
      #                               sep = ""),
      #                 bits.per.sample = 8L, compression = "none",
      #                 reduce = TRUE)
      Image_nuclei_numbers_proteins <- EBImage::Image(data = Image_nuclei_numbers_proteins, colormode = "Color")
      EBImage::writeImage(x = Image_nuclei_numbers_proteins,
                          files = paste(output_dir, image_name_wo_czi, "_nuclei_numbers_proteins1_nuc.tif", sep = ""),
                          type = "tiff",
                          bits.per.sample = 8)
    }

    # Images with marked nuclei and borders around the third protein in
    # cell bodies (proteins in cytosol)
    if(add_scale_bar){
      Image_cytosol_numbers_proteins <- addScaleBar(
        image = Image_cytosol_numbers_proteins,
        length_per_pixel = length_per_pixel_x_in_um)
    }
    if(!is.null(protein_in_cytosol_color) & protein_in_cytosol_color != "none"){
      # tiff::writeTIFF(what = Image_cytosol_numbers_proteins,
      #                 where = paste(output_dir,
      #                               image_name_wo_czi,
      #                               "_nuclei_numbers_proteins2_cell.tif",
      #                               sep = ""),
      #                 bits.per.sample = 8L, compression = "none",
      #                 reduce = TRUE)
      Image_cytosol_numbers_proteins <- EBImage::Image(data = Image_cytosol_numbers_proteins, colormode = "Color")
      EBImage::writeImage(x = Image_cytosol_numbers_proteins,
                          files = paste(output_dir, image_name_wo_czi, "_nuclei_numbers_proteins2_cell.tif", sep = ""),
                          type = "tiff",
                          bits.per.sample = 8)
    }

    # Images with left out nuclei and only with positions of nuclei
    if(add_scale_bar){
      Image_nucleus_part <- addScaleBar(
        image = Image_nucleus_part,
        length_per_pixel = length_per_pixel_x_in_um)
      Image_non_nucleus_part <- addScaleBar(
        image = Image_non_nucleus_part,
        length_per_pixel = length_per_pixel_x_in_um)
    }
    # tiff::writeTIFF(what = Image_nucleus_part,
    #                 where = paste(output_dir,
    #                               image_name_wo_czi,
    #                               "_nucleus_part.tif",
    #                               sep = ""),
    #                 bits.per.sample = 8L, compression = "none",
    #                 reduce = TRUE)

    # (Image_nucleus_part contains also nuclei (parts) that have been deleted
    # for nmask because they are too small or at the borders of the image.)
    Image_nucleus_part <- EBImage::Image(data = Image_nucleus_part, colormode = "Color")
    EBImage::writeImage(x = Image_nucleus_part,
                        files = paste(output_dir, image_name_wo_czi, "_nucleus_part.tif", sep = ""),
                        type = "tiff",
                        bits.per.sample = 8)

    # tiff::writeTIFF(what = Image_non_nucleus_part,
    #                 where = paste(output_dir,
    #                               image_name_wo_czi,
    #                               "_nucleus_left_out.tif",
    #                               sep = ""),
    #                 bits.per.sample = 8L, compression = "none",
    #                 reduce = TRUE)
    Image_non_nucleus_part <- EBImage::Image(data = Image_non_nucleus_part, colormode = "Color")
    EBImage::writeImage(x = Image_non_nucleus_part,
                        files = paste(output_dir, image_name_wo_czi, "_nucleus_left_out.tif", sep = ""),
                        type = "tiff",
                        bits.per.sample = 8)

    # Normalized images with left out clusters and only with positions of clusters
    if(protein_in_membrane_color != "none" & apotome_section){
      if(add_scale_bar){
        Image_clusters_part <- addScaleBar(
          image = Image_clusters_part,
          length_per_pixel = length_per_pixel_x_in_um)
        Image_non_clusters_part <- addScaleBar(
          image = Image_non_clusters_part,
          length_per_pixel = length_per_pixel_x_in_um)
      }
      # tiff::writeTIFF(what = Image_clusters_part,
      #                 where = paste(output_dir,
      #                               image_name_wo_czi,
      #                               "_normalized_clusters_part.tif",
      #                               sep = ""),
      #                 bits.per.sample = 8L, compression = "none",
      #                 reduce = TRUE)
      Image_clusters_part <- EBImage::Image(data = Image_clusters_part, colormode = "Color")
      EBImage::writeImage(x = Image_clusters_part,
                          files = paste(output_dir, image_name_wo_czi, "_normalized_clusters_part.tif", sep = ""),
                          type = "tiff",
                          bits.per.sample = 8)

      # tiff::writeTIFF(what = Image_non_clusters_part,
      #                 where = paste(output_dir,
      #                               image_name_wo_czi,
      #                               "_normalized_clusters_left_out.tif",
      #                               sep = ""),
      #                 bits.per.sample = 8L, compression = "none",
      #                 reduce = TRUE)
      Image_non_clusters_part <- EBImage::Image(data = Image_non_clusters_part, colormode = "Color")
      EBImage::writeImage(x = Image_non_clusters_part,
                          files = paste(output_dir, image_name_wo_czi, "_normalized_clusters_left_out.tif", sep = ""),
                          type = "tiff",
                          bits.per.sample = 8)
    }


    # Remove all variables but the ones used before the for loop

    list_of_variables <- ls()
    keep_variables <- c("df_results", "zis",
                        "bit_depth", "file_names", "image_format",
                        "input_dir",
                        "apotome",
                        "apotome_section",
                        "nucleus_color",
                        "protein_in_nuc_color",
                        "protein_in_cytosol_color",
                        "protein_in_membrane_color",
                        "number_size_factor",
                        "bit_depth",
                        "number_of_pixels_at_border_to_disregard",
                        "add_scale_bar",
                        "thresh_w_h_nuc",
                        "thresh_offset",
                        "thresh_offset_protein_in_nucleus",
                        "thresh_offset_protein_in_cytosol",
                        "blur_sigma",
                        "use_histogram_equalized",
                        "metadata_file",
                        "normalize_nuclei_layer",
                        "magnification_objective",
                        "output_dir",
                        "number_of_images",
                        ".old.options")

    remove_variables <- list_of_variables[
      !(list_of_variables %in% keep_variables)]

    rm(list = remove_variables)
    rm(remove_variables)


  } # end of the routine for every image in the directory
  rm(i)

  if(!is.null(df_results)){
    utils::write.csv(df_results,
                     file = paste(output_dir, "image_analysis_summary_en.csv", sep=""), row.names = FALSE)

    utils::write.csv2(df_results,
                      file = paste(output_dir, "image_analysis_summary_de.csv", sep=""), row.names = FALSE)
  }

  return(df_results)

}

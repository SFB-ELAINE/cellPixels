#' @title addNumberToImage
#' @description Adds a number to the image layer
#' @details Adds an integer to a x-y-3(rgb)-representation of an image.
#' @aliases addnumbertoimage addNumbertoimage addNumberToimage
#' addnumberToImage addnumbertoImage
#' @author Kai Budde
#' @export addNumberToImage
#' @param image An three-dimensional array of numbers between 0 and 1
#' @param number A number (integer to be drawn/copied)
#' @param pos_x A number (x-value for starting pixel (upper left corner))
#' @param pos_y A number (y-value for starting pixel (upper left corner))
#' @param number_size_factor A number (factor for resizing the number)
#' @param number_color A string (color of thr number to be drawn)

addNumberToImage <- function(image = NULL,
                             number = NULL,
                             pos_x = NULL,
                             pos_y = NULL,
                             number_size_factor = NULL,
                             number_color = "red"){

  # Default values for missing arguments -----------------------------------
  if(is.null(image)){
    print(paste("Please call the function with an image.", sep=""))
    return()
  }
  if(is.null(number)){
    print(paste("Please call the function with a number", sep=""))
    return()
  }

  if(is.null(pos_x)){
    pos_x <- 1
  }
  if(is.null(pos_y)){
    pos_y <- 1
  }

  # Read in correct integers -----------------------------------------------
  list_of_digits <- unlist(strsplit(x = toString(number), split = ""))
  number_of_digits <- length(list_of_digits)

  # Add the digits to the image --------------------------------------------
  image_with_numbers <- image

  # Go through all digit images and save the matrices in a list ---
  for(i in 1:number_of_digits){
    digit <- list_of_digits[i]

    digit_path <- paste(digit, ".tiff", sep="")
    digit_path <- system.file("digits", digit_path, package = "cellPixels")

    # digit_image2 <- tiff::readTIFF(source = digit_path, convert = TRUE,
    #                         info = FALSE)
    digit_image <- EBImage::readImage(files = digit_path,type = "tiff")

    # Only keep the black layer
    digit_image <- digit_image[,,4]
    if(i == 1){
      digit_images <- list()

      # Resizing factor depending on size of first digit
      if(is.null(number_size_factor)){
        min_repeating_number <-
          min(dim(image_with_numbers)[1]/dim(digit_image)[1],
              dim(image_with_numbers)[2]/dim(digit_image)[2])

        # It should be possible to have every number at least thirty times
        # in every direction
        number_size_factor <- min_repeating_number / 30
      }


      digits_size_x <- 0
      digits_size_y <- 0
    }

    # Resize image
    digit_image <- resizeImage(image = digit_image,
                               resize_factor = number_size_factor)
    digit_images[[i]] <- digit_image

    # Save dimensions of all digits
    # digits_size_x <- max(digits_size_x, dim(digit_image)[1])
    # digits_size_y <- digits_size_y + dim(digit_image)[2]
    digits_size_x <- digits_size_x + dim(digit_image)[1]
    digits_size_y <- max(digits_size_y, dim(digit_image)[2])
  }

  # Choose layer of image where the number should be added to
  if(tolower(number_color) == "red"){
    image_layer <- 1
  }
  if(tolower(number_color) == "green"){
    image_layer <- 2
  }
  if(tolower(number_color) == "blue"){
    image_layer <- 3
  }
  if(tolower(number_color) == "white"){
    image_layer <- 1:dim(image)[3]
  }

  # Adapt the starting position so the digits will be completely seen ---

  # Right side
  # if( (pos_y + digits_size_y) > dim(image_with_numbers)[2] ){
  #   pos_y <- dim(image_with_numbers)[2] - digits_size_y
  # }
  if( (pos_x + digits_size_x) >= dim(image_with_numbers)[1]){
    pos_x <- dim(image_with_numbers)[1] - digits_size_x
  }

  # Lower side
  # if( (pos_x + digits_size_x) > dim(image_with_numbers)[1]){
  #   pos_x <- dim(image_with_numbers)[1] - digits_size_x
  # }
  if( (pos_y + digits_size_y) >= dim(image_with_numbers)[2] ){
    pos_y <- dim(image_with_numbers)[2] - digits_size_y
  }


  # Add all digit images to bigger one ---
  for(i in 1:number_of_digits){
    digit_image_size_x <- dim(digit_images[[i]])[1]
    digit_image_size_y <- dim(digit_images[[i]])[2]


    # for(row in 1:digit_image_size_x){
    #   for(col in 1:digit_image_size_y){
    #     image_with_numbers[pos_x + row, pos_y + col, image_layer] <-
    #       min(digit_images[[i]][row, col] +
    #             image_with_numbers[pos_x + row, pos_y + col, image_layer], 1)
    #   }
    # }
    for(col in 1:digit_image_size_x){
      for(row in 1:digit_image_size_y){
        image_with_numbers[pos_x + col, pos_y + row, image_layer] <-
          min(digit_images[[i]][col, row] +
                image_with_numbers[pos_x + col, pos_y + row, image_layer], 1)
      }
    }

    # Adapt starting position for the next digit
    #pos_y <- pos_y + digit_image_size_y
    pos_x <- pos_x + digit_image_size_x
  }

  return(image_with_numbers)
}

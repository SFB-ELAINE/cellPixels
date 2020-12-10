#' @title resizeImage
#' @description Resizes an image
#' @details Resizes an image with a given scaling factor
#' @aliases resizeimage
#' @author Kai Budde
#' @export resizeImage
#' @param image An one to three-dimensional array of numbers between 0 and 1
#' @param number_size_factor A number (factor for resizing the number)

resizeImage <- function(image = NULL,
                        number_size_factor = 1){

  # Default values for missing arguments -----------------------------------
  if(is.null(image)){
    print(paste("Please call the function with an image.", sep=""))
    return()
  }
  if(is.null(number_size_factor)){
    print(paste("Please call the function with a number", sep=""))
    return()
  }

  if(number_size_factor < 0){
    print(paste("The factor may only be a number > 0", sep=""))
    return()
  }

  # Leave image as is if factor == 1 ---------------------------------------
  if(number_size_factor == 1){
    return(image)
  }

  # Shrink image if number_size_factor </> 1 ---------------------------------
  if(number_size_factor > 0){

    new_height <-  round(number_size_factor * dim(image)[1])
    new_width  <-  round(number_size_factor * dim(image)[2])

    new_img = apply(image, 2, function(y){return (stats::spline(y, n = new_height)$y)})
    new_img = t(apply(new_img, 1, function(y){return (stats::spline(y, n = new_width)$y)}))

    new_img[new_img < 0] = 0
    new_img = round(new_img)

    return(new_img)
  }

  print("Something went wrong.")
  return(0)
}

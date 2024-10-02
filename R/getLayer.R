#' @title getLayer
#' @description Get a specific layer of the image
#' @details By using an x-y-3(rgb)-representation of an image, you can
#' extract one layer.
#' @aliases getlayer
#' @author Kai Budde-Sagert
#' @export getLayer
#' @param image An three-dimensional array of numbers between 0 and 1
#' @param layer A character (color of the layer)

getLayer <- function(image = NULL,
                     layer = NULL){

  # Default values for missing arguments -----------------------------------
  if(is.null(image)){
    print(paste("Please call the function with an image.", sep=""))
    return()
  }

  if(missing(layer)){
    layer <- "red"
  }

  layer <- tolower(layer)

  # Extract the layer of the cilia -----------------------------------------

  if(layer == "red" | layer == "r"){
    image_layer <- image[,,1]
  }else if(layer == "green" | layer == "g"){
    image_layer <- image[,,2]
  }else if(layer == "blue" | layer == "b"){
    image_layer <- image[,,3]
  }else{
    print(paste("Please enter a color (red, green or blue).",
                sep=""))
    return()
  }

  return(image_layer)
}


#' Creating a Neighborhood Matrix for Areal Data from a Shapefile
#' @description Takes a path to a shape file and creates a neighborhood matrix.  This neighborhood
#' matrix can be used with the objective ICAR model \insertCite{Keefe_2018}{ref.ICAR}.
#'
#' @author Erica M. Porter, Matthew J. Keefe, Christopher T. Franck, and Marco A.R. Ferreira
#'
#' @param shape.file File path to a shapefile.
#' @importFrom sf st_read
#' @importFrom spdep poly2nb nb2mat
#' @return A list containing a neighborhood matrix and the SpatialPolygonsDataFrame object
#' corresponding to the shape file.
#'     \item{H}{A neighborhood matrix.}
#'     \item{map}{SpatialPolygonsDataFrame object from the provided shapefile.}
#' @export
#' @examples
#' #### Load extra libraries
#' library(sp)
#' library(sf)
#'
#' ### Read in a shapefile of the contiguous U.S. from package data
#' system.path <- system.file("extdata", "us.shape48.shp", package = "ref.ICAR", mustWork = TRUE)
#' shp.layer <- gsub('.shp','',basename(system.path))
#' shp.path <- dirname(system.path)
#' us.shape48 <- st_read(dsn = path.expand(shp.path), layer = shp.layer)
#'
#' shp.data <- shape.H(system.path)
#' names(shp.data)
shape.H <- function(shape.file) {

  shp.layer <- gsub('.shp','',basename(shape.file))
  #shp.path <- gsub(basename(shape.file),'',shape.file)
  shp.path <- dirname(shape.file)
  map <- st_read(dsn = path.expand(shp.path), layer = shp.layer)

  map.nb <- poly2nb(map)
  W <- nb2mat(map.nb, style = "B")
  Dmat <- diag(apply(W,1,sum))
  lambda <- 0

  H <- Dmat - W
  row.names(H) <- NULL
  return(list(H=H,map=map))
}

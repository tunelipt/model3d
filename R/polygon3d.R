#' Create a new simple 3D polygon.
#'
#' Creates a new simple 3D polygon from the coordinates. The polygon is stored as a matrix where the first row are the x coordinates of the vertices and the second row are the y coordinates anf the third row are the z coordinates.
#'
#' @param x Vector with x coordinates of the vertices
#' @param y Vector with y coordinates of the vertices
#' @param z Vector with z coordinates of the vertices
#' @return A matrix with class \code{polygon3d} where each column is a vertex of the polygon
#' @seealso \code{\link{newPolygon}}, \code{\link{polygonArea}}, \code{\link{polygonCentroid}}
#' @export
#' @examples
#' p <- newPolygon3d(c(0,1,2,1,2,0), c(0,0,1,2,3,3), c(0,0,0,0,0,0))
#' print(p)
newPolygon3d <- function(x, y, z){
  p <- rbind(x=x, y=y, z=z)
  class(p) <- c('polygon3d', 'polygon')
  return(p)
}

#' Calculates the sum of the cross product of all points of a polygon.
#'
#' Helper function that calculates the sum of the cross products between
#' sequential vector formed between the ith point and the first.
#'
#' @export
polyCross <- function(p){
  np <- dim(p)[2]

  p1 <- p[,1]

  area <- c(0,0,0)

  for (i in 3:np)
    area <- area + crossProduct(p[,i-1] - p1, p[,i]-p1)

  return(as.double(area))
}

  

#' Calculates the area of a polygon.
#'
#' A method that computes the area of simple convex 3D polygons of class polygon3d.
#'
#' @seealso \code{\link{polygonArea}}
#' @export
#' @examples
#' p <- newPolygon3d(c(0,1,2,1,2,0), c(0,0,1,2,3,3), c(0,0,0,0,0,0))
#' print(polygonArea(p))
polygonArea.polygon3d <- function(p){

  area <- polyCross(p)

  return(vnorm(area)/2)
  
}


#' Returns the normal to a polygon.
#'
#' Assumes the polygon is convex.
#' @export
polygonNormal <- function(p) UseMethod("polygonNormal")

#' Returns the normal to a polygon.
#'
#' For 2D polygons, just return +1 or -1, depending on orientation.
#' @export
polygonNormal.polygon <- function(p) polygonOrient(p)*c(0,0,1)


#' Returns the normal to a polygon.
#'
#' Method for polygon3d objects. Assumes the polygon is convex.
#' @export
polygonNormal.polygon3d <- function(p){
  

  area <- polyCross(p)

  return(as.double(area / vnorm(area)))
  
}         


#' Centroid of a 3D polygon.
#'
#' Method that computes centroid of a polygons of class polygon3d.
#'
#' @param p Polygon
#' @return Coordinates of the centroid
#' @export
polygonCentroid.polygon3d <- function(p){


  np <- dim(p)[2]

  p1 <- p[,1]

  area <- matrix(0, nr=3, nc=np-2)
  r <- matrix(0, nr=3, nc=np-2)

  for (i in 3:np)
    area[,i-2] <- 0.5 * crossProduct(p[,i-1] - p1, p[,i]-p1)

  A <- rowSums(area)
  AA <- vnorm(A)
  nn <- A / AA


  for (i in 3:np)
    r[,i-2] <- (p[,1] + p[,i-1] + p[,i]) / 3 * sum(area[,i-2]*nn)

  return( rowSums(r) / AA)
}

  
  

#' Method for 3D polygon orientation.
#'
#' The function returns the normal to the polygon
#' according to right hand rule orientation.
#' If the area is zero, return the null vector c(0,0,0).
#' @export
polygonOrient.polygon3d <- function(p){

  area <- polyCross(p)
  
  v <- vnorm(area)

  if (abs(v) > 0)
    return(as.double(area/v))
  return(c(0,0,0))
}



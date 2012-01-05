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
  n <- polygonNormal.polygon3d(p)
  
  np <- dim(p)[2]
  p <- cbind(p, p[,1])

  accu <- double(3)

  for (i in 1:np)
    accu <- accu + crossProduct(p[,i], p[,i+1])

  A <- abs(dotProduct(n, accu)) / 2
  return(A)
}


#' Returns the normal to a polygon.
#'
#' Assumes the polygon is convex.
#' @export
polygonNormal <- function(p) UseMethod("polygonNormal")

#' Returns the normal to a polygon.
#'
#' For 2D polygons, just return z.
#' @export
polygonNormal.default <- function(p) c(0, 0, 1)


#' Returns the normal to a polygon.
#'
#' Method for polygon3d objects. Assumes the polygon is convex.
#' @export
polygonNormal.polygon3d <- function(p){
  
  np <- dim(p)[2]
  if (np<3) return(NULL)

  cp <- crossProduct( p[,2]-p[,1], p[,3]-p[,1] )
  return(as.double(cp / vnorm(cp)))
  
}         




require(gpclib)
require(deldir)

#' Create a new simple 2D polygon.
#'
#' Creates a new simple 2D polygon from the coordinates. The polygon is stored as a matrix where the first row are the x coordinates of the vertices and the second row are the y coordinates.
#'
#' @param x Vector with x coordinates of the vertices
#' @param y Vector with y coordinates of the vertices
#' @return A matrix with class \code{polygon} where each column is a vertex of the polygon
#' @seealso \code{\link{newPolygon3d}}, \code{\link{polygonArea}}, \code{\link{polygonCentroid}}, \code{\link{polygonIntersect}}
#' @export
#' @examples
#' p <- newPolygon(c(0,1,2,1,2,0), c(0,0,1,2,3,3))
#' print(p)
newPolygon <- function(x, y){
  p <- rbind(x=x, y=y)
  class(p) <- 'polygon'
  return(p)
}

#' Converts a polygon to a gpclib polygon.
#'
#' The R package gpclib provides a few algorithms that are useful but it has an internal (very general) representation of polygons. This function converts a simple polygon of class polygon to a gpclib polygon.
#'
#' @param p A polygon.
#' @return A gpclib polygon.
#' @seealso \code{\link{newPolygon}}, \code{\link{gpc2poly}}
#' @export
#' @examples
#' p <- newPolygon(c(0,1,2,1,2,0), c(0,0,1,2,3,3))
#' g <- poly2gpc(p)
#' plot(g)
poly2gpc <- function(p) as(cbind(p[1,], p[2,]), 'gpc.poly')

#' Converts a gpclib polygon to a simple polygon.
#'
#' This function Assumes that the polygon is simple and has no holes.
#'
#' @param p A gpclib polygon
#' @return A polygon of class polygon
#' @seealso \code{\link{newPolygon}}, \code{\link{poly2gpc}}
#' @export
#' @examples
#' p <- newPolygon(c(0,1,2,1,2,0), c(0,0,1,2,3,3))
#' g <- poly2gpc(p)
#' print(gpc2poly(g))
gpc2poly <- function(gpc){p <- get.pts(gpc)[[1]]; return(newPolygon(p$x, p$y))}


#' Returns the coordinates of a simple polygon.
#'
#' Returns a matrix where each row is a vertex of the polygon. Only works with simples polygons. It is actually a generic S3 function.
#'
#' @param p Polygon
#' @return Matrix with coordinates of the vertices
#' @export
getCoords <- function(p) UseMethod("getCoords")

#' Method of getCoords for polygon.
#' @export
getCoords.polygon <- function(p){
  x <- t(p)
  class(x) <- 'matrix'
  return(x)
}

#' Method of getCoords for gpc.poly.
#' @export
getCoords.gpc.poly <- function(p){
  x <- get.pts(p)[[1]]

  return(cbind(x=x$x, x=x$y))
}
         
#' Calculates the area of a polygon.
#'
#' A generic S3 function that computes the area of polygons given the coordinates of its vertices arranged as a matrix where each column is a vertex numbered in counterclockwise manner.
#'
#' @param p Matrix containing the vertices of the polygon.
#' @return Area of the polygon.
#' @seealso \code{\link{crossProduct}}, \code{\link{poly3dArea}}
#' @export
#' @examples
#' poly <- newPolygon(c(0, 1, 1, 0), c(0, 0, 1, 1))
#' print(polygonArea(poly))
polygonArea <- function(p) UseMethod("polygonArea")

#' Calculates the area of a polygon.
#'
#' A method that computes the area of simple polygons of class polygon.
#'
#' @seealso \code{\link{polygonArea}}, \code{\link{polygonArea.gpc.poly}}
#' @export
#' @examples
#' poly <- newPolygon(c(0, 1, 1, 0), c(0, 0, 1, 1))
#' print(polygonArea(poly))
polygonArea.polygon <- function(p){
  n <- dim(p)[2]
  x <- c(p[1,], p[1,1])
  y <- c(p[2,], p[2,1])
  area <- (x[1:n]*y[2:(n+1)] - x[2:(n+1)]*y[1:n])
  area <- 0.5*(sum(area))
  return(area)
}


#' Calculates the area of a polygon.
#'
#' A method that computes the area of polygons of class gpc.poly.
#'
#' @seealso \code{\link{polygonArea}}, \code{\link{polygonArea.polygon}}
#' @export
#' @examples
#' p <- newPolygon(c(0,1,2,1,2,0), c(0,0,1,2,3,3))
#' g <- poly2gpc(p)
#' print(polygonArea(g)
polygonArea.gpc.poly <- function(p) area.poly(p)




#' Centroid of a polygon.
#'
#' Generic S3 function that computes the centroid of a polygon.
#'
#' @param p Polygon
#' @return Coordinates of the centroid
#' @export
polygonCentroid <- function(p) UseMethod("polygonCentroid")


#' Centroid of a polygon.
#'
#' Method that computes centroid of a polygons of class polygon.
#'
#' @param p Polygon
#' @return Coordinates of the centroid
#' @export
polygonCentroid.polygon <- function(p){
  A <- polygonArea(p)
  x <- p[1,]
  y <- p[2,]
  n <- length(x)
  x <- c(x, x[1])
  y <- c(y, y[1])
  
  cx <- (x[1:n]+x[2:(n+1)]) * (x[1:n]*y[2:(n+1)] - x[2:(n+1)]*y[1:n])
  cx <- 1/(6*A) * sum(cx)
  cy <- (y[1:n]+y[2:(n+1)]) * (x[1:n]*y[2:(n+1)] - x[2:(n+1)]*y[1:n])
  cy <- 1/(6*A) * sum(cy)

  return(c(cx, cy))
}
  


#' Polygon Intersection.
#'
#' S3 generic function that compute the intersection of two polygons. Remember that the intersection of two polygons can be 0, 1 or more polygons.
#'
#' @param p1 First polygon
#' @param p2 Second polygon
#' @return The intersection of the polygons.
#' @export
polygonIntersect <- function(p1, p2) UseMethod('polygonIntersect')

#' Polygon Intersection.
#'
#' Method that computes  intersection of two gpc.poly polygons. Actually calls the approppriate function of gpclib to compute the intersection.
#' @export
polygonIntersect.gpc.poly <- function(p1, p2){

  return(gpclib::intersect(p1, p2))
}
#' Polygon Intersection.
#'
#' Method that computes  intersection of two polygons of class polygon. It converts the polygons to gpc.poly polygons and then uses gpclib to compute the intersection. The result is then converted to a list of polygons of class polygon.
#' @export
polygonIntersect.polygon <- function(p1, p2){

  g1 <- poly2gpc(p1)
  g2 <- poly2gpc(p2)

  g <- gpclib::intersect(g1, g2)

  p <- lapply(get.pts(g), function(x) newPolygon(x$x, x$y))
  return(p)
}



  

#' Centroid of a polygon.
#'
#' Method that computes centroid of a polygons of class goc.poly.
#'
#' @param p Polygon
#' @return Coordinates of the centroid
#' @export
polygonCentroid.gpc.poly <- function(p){
  pts <- get.pts(p)
  np <- length(pts)
  pp <- lapply(pts, function(p) newPolygon(p$x, p$y))
  
  A <- sapply(pp, function(p) polygonArea.polygon(p))
  C <- sapply(pp, function(p) polygonCentroid.polygon(p))
  #hole <- sapply(pts, function(p) p$hole)
  #fhole <- double(np)
  #fhole[hole] <- -1
  #fhloe[!hole] <- 1
  
  xc <- sum(C[1,] * A) / sum(A)
  yc <- sum(C[2,] * A) / sum(A)

  return(c(xc, yc))
}
  
  



  
  


               
     


  
  


#' Computes the voronoi tesselation of a set of 2D points.
#'
#' Uses the package \code{\link{deldir}} to compute the voronoi tesselation of a set of points.
#'
#' @param x X coordinates of the points
#' @param y Y coordinates of the points
#' @param xb X limits of outer box of the tesselation
#' @param yb Y limits of outer box of the tesselation
#' @return A list of polygons corresponding to the voronoi regions of each node.
#' @export
#' @examples
#' x <- runif(10)
#' y <- runif(10)
#' vor <- voronoi2d(x, y, c(0, 1), c(0, 1))
#' plot(x, y, ty='p', xlim=c(-.1, 1.1), ylim=c(-.1, 1.1))
#' for (p in vor) polygon(v[1,], v[2,])
voronoi2d <- function(x, y, xb=NULL, yb=NULL){
  if (length(x)==1){
    if (is.null(xb)) xb <- c(x-10, x+10)
    if (is.null(yb)) yb <- c(y-10, y+10)

    xx <- c(xb[1], xb[2], xb[2], xb[1])
    yy <- c(yb[1], yb[1], yb[2], yb[2])

    return(list(newPolygon(xx, yy)))
  }

  if (is.null(xb)){
    rx <- range(x)
    dx <-rx[2] - rx[1]
    xb <- c(rx[1] - 1.2*dx, rx[2] + 1.2*dx)
  }

  if (is.null(yb)){
    ry <- range(y)
    dy <-ry[2] - ry[1]
    yb <- c(ry[1] - 1.2*dy, ry[2] + 1.2*dy)
  }


    
  tri <- deldir(x,y, rw=c(xb[1], xb[2], yb[1], yb[2]))
  vor <- tile.list(tri)
  plst <- lapply(vor, function(v) newPolygon(v$x, v$y))
  
  return(plst)

}


  

#' Is the point inside the the polygon?
#'
#' Determines if a point is inside a simple polygon.
#'
#' @param vx x coordinates of the vertices of the polygon.
#' @param vy y coordinates of the vertices of the polygon.
#' @param x x coordinate of the point to be tested.
#' @param y y coordinate of the point to be tested.
#' @return TRUE if the point is inside the polygon.
#' @export
#' @examples
#' vx <- c(0, 2, 2, 0)
#' vy <- c(0, 0, 2, 2)
#' print(pnpoly(vx, vy, 1, 1))
#' print(pnpoly(vx, vy, 3, 1))
pnpoly <- function(vx, vy, x, y){
  nv <- length(vx)
  test <- FALSE

  j <- nv
  for (i in 1:nv){
    if ( ((vy[i]>y) != (vy[j]>y)) &&
         ( x <= (vx[j]-vx[i]) * (y - vy[i]) / (vy[j] - vy[i]) + vx[i]))
      test <- !test
    j <- i
  }
  return(test)
  
}

#' Is the point in the boundary of a polygon?
#'
#' Determines if a point is on the boundary of a polygon.
#'
#' @param vx x coordinates of the vertices of the polygon.
#' @param vy y coordinates of the vertices of the polygon.
#' @param x x coordinate of the point to be tested.
#' @param y y coordinate of the point to be tested.
#' @param eps Admissible error.
#' @return TRUE if the point is on the segments.
#' @export
#' @examples
#' vx <- c(0, 2, 2, 0)
#' vy <- c(0, 0, 2, 2)
#' print(pnpolybnd(vx, vy, 1, 1))
#' print(pnpolybnd(vx, vy, 2, 1))
pnpolybnd <- function(vx, vy, x, y, eps0=1e-8){
  vx <- c(vx, vx[1])
  vy <- c(vy, vy[1])
  
  n <- length(vx)
  nm1 <- n-1

 
 for (i in 2:n){
   im1 <- i-1
   dx <- vx[i] - vx[im1]
   dy <- vy[i] - vy[im1]
   alfx=NA
   alfy=NA
   
   if (abs(dy) > eps0)
     alfy <- (y - vy[im1]) / (vy[i] - vy[im1])
   if (abs(dx) > eps0)
     alfx <- (x - vx[im1]) / (vx[i] - vx[im1])
   
   if (is.na(alfy)){
     if(abs(y-vy[i])<eps0 && alfx >= 0 && alfx <= 1)
       return(TRUE)
   }else if (is.na(alfx)){
     if (abs(x-vx[i])<eps0 && alfy >= 0 && alfy <= 1)
       return(TRUE)
   }else{
     if (abs(alfx-alfy) < eps0) return(TRUE)
   }
 }

 return(FALSE)
}



#' Is the point inside or on the boundary of a polygon?
#'
#' Determines if a point is inside or on the boundary of a simple polygon.
#' This function simply calls \code{\link{pnpoly}} and \code{\link{pnpolybnd}}.
#'
#' @param vx x coordinates of the vertices of the polygon.
#' @param vy y coordinates of the vertices of the polygon.
#' @param x x coordinate of the point to be tested.
#' @param y y coordinate of the point to be tested.
#' @return TRUE if the point is inside the polygon.
#' @export
#' @examples
#' vx <- c(0, 2, 2, 0)
#' vy <- c(0, 0, 2, 2)
#' print(pnpoly(vx, vy, 1, 1))
#' print(pnpoly(vx, vy, 3, 1))
#' print(pnpoly(vx, vy, 2, 1))
ponpoly <- function(vx, vy, x, y, eps0=1e-8){

  return(pnpoly(vx, vy, x, y) || pnpolybnd(vx, vy, x, y, eps0))
}

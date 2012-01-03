require(gpclib)
require(deldir)


newPolygon <- function(x, y){
  p <- rbind(x=x, y=y)
  class(p) <- 'polygon'
  return(p)
}


poly2gpc <- function(p) as(cbind(p[1,], p[2,]), 'gpc.poly')

gpc2poly <- function(gpc){p <- get.pts(gpc)[[1]]; return(newPolygon(p$x, p$y))}


getCoords <- function(p) UseMethod("getCoords")

getCoords.polygon <- function(p){
  x <- t(p)
  class(x) <- 'matrix'
  return(x)
}

getCoords.gpc.poly <- function(p){
  x <- get.pts(p)[[1]]

  return(cbind(x=x$x, x=x$y))
}
         
polygonArea <- function(p) UseMethod("polygonArea")

polygonArea.polygon <- function(p){
  n <- dim(p)[2]
  x <- c(p[1,], p[1,1])
  y <- c(p[2,], p[2,1])
  area <- (x[1:n]*y[2:(n+1)] - x[2:(n+1)]*y[1:n])
  area <- 0.5*(sum(area))
  return(area)
}


polygonArea.gpc.poly <- function(p) area.poly(p)





polygonCentroid <- function(p) UseMethod("polygonCentroid")

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
  


splitWith <- function(vor, hull){

  return(lapply(vor, function(v) intersect(v, hull)))
}


polygonIntersect <- function(p1, p2) UseMethod('polygonIntersect')

polygonIntersect.gpc.poly <- function(p1, p2){

  return(intersect(p1, p2))
}
polygonIntersect.polygon <- function(p1, p2){

  g1 <- poly2gpc(p1)
  g2 <- poly2gpc(p2)

  g <- intersect(g1, g2)

  p <- lapply(get.pts(g), function(x) newPolygon(x$x, x$y))
  return(p)
}



  

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
  
  



  
  


               
     


  
  



voronoi2d <- function(x, y, xb=NULL, yb=NULL){
  if (length(x)==1){
    if (is.null(xb)) xb <- c(x-10, x+10)
    if (is.null(yb)) yb <- c(y-10, y+10)

    xx <- c(xb[1], xb[2], xb[2], xb[1])
    yy <- c(yb[1], yb[1], yb[2], yb[2])

    return(list(as(cbind(xx, yy), 'gpc.poly')))
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


  


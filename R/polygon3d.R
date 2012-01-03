newPolygon3d <- function(x, y, z){
  p <- rbind(x=x, y=y, z=z)
  class(p) <- c('polygon3d', 'polygon')
  return(p)
}


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


polygonNormal <- function(p) UseMethod("polygonNormal")
polygonNormal.default <- function(p) c(0, 0, 1)


polygonNormal.polygon3d <- function(p){
  
  np <- dim(p)[2]
  if (np<3) return(NULL)

  cp <- crossProduct( p[,2]-p[,1], p[,3]-p[,1] )
  return(as.double(cp / vnorm(cp)))
  
}         




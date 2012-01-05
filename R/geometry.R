                                        # Geometru operations:


#' Cross product of two vectors.
#'
#' Compute the cross product of two 3D vectors.
#'
#' @param u First vector
#' @param v Second vector
#' @return u X v
#' @seealso \code{\link{dotProduct}}
#' @examples
#' ex <- c(1, 0, 0)
#' ey <- c(0, 1, 0)
#' ez <- crossProduct(ex, ey)
crossProduct <- function(u,v){

  x1 <- u[2]*v[3] - u[3]*v[2]
  x2 <- u[3]*v[1] - u[1]*v[3]
  x3 <- u[1]*v[2] - u[2]*v[1]
  return(c(x1, x2, x3))
}


#' Dot product of two vectors.
#'
#' Compute the dot product of two n-dimensional
#' vectors.
#'
#' @param u First vector
#' @param v Second vector
#' @return \eqn{\sum_{i=1}^{n} u_i v_i}
#' @seealso \code{\link{crossProduct}}, \code{\link{vnorm}}
#' @examples
#' print(dotProduct(c(1,0,0), c(0,1,0)))
dotProduct <- function(u,v)
  return(sum(u*v))

#' Euclidean norm of a vector.
#'
#' Compute the norm of a n-dimensional vector.
#'
#' @param u Vector
#' @return \eqn{\sqrt{\sum_{i=1}^{n} u_i v_i}}
vnorm <- function(u)
  return(sqrt(dotProduct(u,u)))

#' Compute the normal to a triangle.
#'
#' Compute the normal to the plane formed by the edges of a triangle.
#'
#' @param x Matrix 3x3 with each column as a vertex.
#' @param invert A boolean that specifies whether the outer normal is wanted.
#' @return A vector containing the 3 components of the normal.
#' @seealso \code{\link{crossProduct}}
#' @examples
#' tri <- matrix(c(0,0,0, 1,0,0, 1,1,0), 3, 3)
#' print(triNormal(tri))
triNormal <- function(x, invert=FALSE){
  u <- x[,2] - x[,1]
  v <- x[,3] - x[,1]

  b <- crossProduct(u,v)
  s <- ifelse(invert, -1, 1)
  return(s*b / vnorm(b))
}

#' Normal of a 3D polygon.
#'
#' Computes the normal of a simple polygon.
#'
#' @param p Polygon where each column are the coordinates of a vertex
#' @param eps Relative error acceptable when determining if 3 vertices are colinear
#' @return A vector containing the normal to the polygon when the vertices are assumed to be numbered in counterclockwise manner.
#' @seealso \code{\link{crossProduct}}, \code{\link{triNormal}}
#' @examples
#' poly <- matrix(c(0,0,0, 1,0,0, 1,1,0, 0,1,0), 3, 4)
#' print(poly3dNorm(poly))
poly3dNorm <- function(p, eps=1e-8){
  np <- dim(p)[2]
  if (np<3) return(NULL)
  ll <- c(max(p[1,]) - min(p[1,]), max(p[2,]) - min(p[2,]), max(p[3,])-min(p[3,]))
  l <- vnorm(ll)

  p1 <- p[,2]-p[,1]

  
  for (i in 3:np){
    cp <- crossProduct(p[,2]-p[,1], p[,i]-p[,1])
    ll <- vnorm(cp)
    if (ll > l*eps)
      return(cp/ll)
  }
  return(c(0,0,0))
  
}

#' Calculates the area of a 3D polygon.
#'
#' Computes the area of a simple 3D polygon given the coordinates of its vertices arranged as a matrix where each column is a vertex numbered in counterclockwise manner.
#'
#' @param p Matrix containing the vertices of the polygon.
#' @param n Normal to the polygon (if it has been previously computed).
#' @return Area of the polygon.
#' @seealso \code{\link{crossProduct}}, \code{\link{poly3dNormal}}
#' @examples
#' poly <- matrix(c(0,0,0, 1,0,0, 1,1,0, 0,1,0), 3, 4)
#' print(poly3dArea(poly))
poly3dArea <- function(p, n=NULL){

  if (is.null(n)) n <- poly3dNorm(p)
  
  np <- dim(p)[2]
  p <- cbind(p, p[,1])

  accu <- double(3)

  for (i in 1:np)
    accu <- accu + crossProduct(p[,i], p[,i+1])

  A <- abs(dotProduct(n, accu)) / 2
  return(A)
    
}


# Projects a point p on a plane given by a normal n and a point p0.
projectPoint <- function(p, n, p0=NULL, eps=1e-6){

  if (is.null(p0)) p0 <- double(3)
  # Verify if p lies on the plane:
  
  ll <- p - p0
  nl <- vnorm(ll)
  lref <- max(nl, 1)

  if (abs(dotProduct(ll, n)) < eps) return(p)
  
  
  A <- matrix(c(1, 0, 0, n[1],
                0, 1, 0, n[2],
                0, 0, 1, n[3],
                n[1], n[2], n[3], 0), 4, 4, byrow=TRUE)
  b <- c(ll, 0)
  x <- solve(A, b)

  return(x[1:3]+p0)
}



generate2dCoords <- function(p){

  x3 <- poly3dNorm(p)
  x1 <- p[,2] - p[,1]
  x1 <- x1 / vnorm(x1)

  x2 <- crossProduct(x3, x1)

  return(cbind(x1, x2, x3))
}



  

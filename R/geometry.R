# Geometru operations:


# Cross product between two vectors
crossProduct <- function(u,v){

  x1 <- u[2]*v[3] - u[3]*v[2]
  x2 <- u[3]*v[1] - u[1]*v[3]
  x3 <- u[1]*v[2] - u[2]*v[1]
  return(c(x1, x2, x3))
}


# Dot product between two vectors
dotProduct <- function(u,v)
  return(sum(u*v))

# Euclidean norm of a vector
vnorm <- function(u)
  return(sqrt(dotProduct(u,u)))

# Normal vector to a triangle
triNormal <- function(x, invert=FALSE){
  u <- x[,2] - x[,1]
  v <- x[,3] - x[,1]

  b <- crossProduct(u,v)
  s <- ifelse(invert, -1, 1)
  return(s*b / vnorm(b))
}

# Determines the normal of a 3D polygon on a plane.
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

# Calculates the area of a 3D polygon.
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



  

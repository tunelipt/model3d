getBase <- function(bz){
  ex <- c(1,0,0)
  ey <- c(0,1,0)
  ez <- c(0,0,1)
  ee <- cbind(x=ex, y=ey, z=ez)
  dd <- sapply(1:3, function(i) dotProduct(bz, ee[,i]))

  
  idx <- which.min(abs(dd))
  py <- ee[,idx]
  by <- py - sum(py*bz) * bz
  by <- by / sqrt(sum(by*by))
  bx <- crossProduct(by, bz)
  base <- cbind(x=bx, y=by, z=bz)
  
  return(base)
}


pointsProject <- function(p, base2d) t(base2d) %*% p




  
polygonProject <- function(p, base2d){
  p2 <- pointsProject(p, base2d)
  class(p2) <- 'polygon'
  return(p2)
}


  
faceProject2d <- function(f, base2d) UseMethod('faceProject2d')

faceProject2d.default <- function(f, base2d) pointsProject(f, base2d)

faceProject2d.polygon3d <- function(f, base2d) polygonProject(f, base2d)

faceProject2d.mesh3d <- function(f, base2d){
  f2 <- apply(f, c(2,3), function(x) pointsProject(x, base2d))
  class(f2) <- 'mesh'
  return(f2)
}

faceProject2d.pmesh3d <- function(f, base2d){
  f2 <- lapply(f, function(p) faceProject2d.polygon3d(p, base2d))
  class(f2) <- 'pmesh'
  return(f2)
}


  
projectionInfo <- function(p, nn, base2d){

  A <- solve(rbind(t(base2d), nn))[,1:2]
  k <- sum(p*nn)
  pinfo <- list(p=p, norm=nn, A=A, base2d=base2d, k=k)
  class(pinfo) <- 'projectionInfo'
  return(pinfo)
  
}



project3dTo2d <- function(p, base2d) UseMethod('project3dTo2d')
project3dTo2d.default <- function(p, base2d) faceProject2d.default(p, base2d)
project3dTo2d.polygon3d <- function(p, base2d) faceProject2d.polygon3d(p, base2d)
project3dTo2d.mesh3d <- function(p, base2d) faceProject2d.mesh3d(p, base2d)
project3dTo2d.pmesh3d <- function(p, base2d) faceProject2d.pmesh3d(p, base2d)


project2dTo3d <- function(p, pinfo) UseMethod('project2dTo3d')
project2dTo3d.default <- function(p, pinfo) pinfo$A %*% c(p, pinfo$k)
project2dTo3d.polygon <- function(p, pinfo){
  p3 <- pinfo$A %*% rbind(p, pinfo$k)
  class(p3) <- 'polygon3d'
  return(p3)
}

project2dTo3d.gpc.poly <- function(p, pinfo) project2dTo3d.polygon(gpc2poly(p), pinfo)



  

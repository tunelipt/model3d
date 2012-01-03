faceProjectionInfo <- function(face, base2d, facor=1){

  nf <- meshSize(face)
  norm <- facor * meshNormal(face)

  pinfo <- lapply(1:nf,
                  function(i){
                    p <- meshGet(face, i)
                    p0 <- p[,1]
                    nn <- norm[,i]
                    return(projectionInfo(p0, nn, base2d))})
  
  return(pinfo)
}

  

surfaceVoronoi <- function(m, pts, basefun=getBase){

  
  x <- pts$x
  y <- pts$y
  z <- pts$z


  mean.norm <- meshMeanNormal(m)

  base2d <- basefun(mean.norm)[,1:2]

  pinfo <- faceProjectionInfo(m, base2d)
  

  pts <- pointsProject(rbind(x,y,z), base)
  npts <- length(x)
  m2d <- project3dTo2d(m, base2d)

  rx <- range(m2d[1,,])
  ry <- range(m2d[2,,])
  dx <- rx[2]-rx[1]
  dy <- ry[2]-ry[1]
  rx[1] <- rx[1]-dx; rx[2] <- rx[2]+dx
  ry[1] <- ry[1]-dy; ry[2] <- ry[2]+dy
  
  vor <- voronoi2d(pts[1,], pts[2,], rx, ry)
  nvor <- length(vor)
  np <- meshSize(m)

  plst <- list()
  for (i in 1:nvor){
    count <- 0
    tmp <- list()
    for (k in 1:np){

      pi <- polygonIntersect(vor[[i]], meshGet(m2d, k))
      if (length(pi) == 0)
        next
      count <- count+1
      pi3d <- project2To3d(pi[[1]], pinfo[[k]])
      tmp[[count]] <- pi3d
    }
    class(tmp) <- 'pmesh3d'
    plst[[i]] <- tmp
  }

  return(plst)
}


    
    

  

  
  
  




intersectHeight <- function(z, h){
  n <- length(z)
  
  dz <- z - h
  x <- dz == 0
  ii <- integer(n)
  ii[x] <- 2
  sig <- dz[1:n]*c(dz[2:n], dz[1])

  ii[sig < 0] <- 1
  return(ii)
}
calcIntersect <- function(z1, z2, z){

  dz <- z2-z1

  if (dz==0) return(NULL)

  return( (z-z1)/dz )
}


  
  
projectIntoHeight <- function(p, hmin, hmax, dnum=2){
  z <- p[dnum,]
  zg <- z >= hmax
  if (all(zg)) return(NULL)
  zs <- z <= hmin
  if (all(zs)) return(NULL)

  # Existe uma intersecção: verificar se está totalmente dentro
  if (all(!zg) && all(!zs)) return(p)
  
  # Sobram 3 casos:
  #
  if (any(zg)){
    n <- length(z)
    idx <- c(1:n,1)
    ii <- intersectHeight(z,hmax)
    ip <- z-hmax <= 0
    pp <- matrix(NA, nr=3, nc=2*n)
    count <- 1
    pt <- 1
    for (pt in 1:n){
      if (ip[pt]){
        pp[,count] <- p[,pt]
        count <- count + 1
      }
      if (ii[pt]==1){
        npt <- idx[pt+1]
        x <- calcIntersect(z[pt], z[npt], hmax)
        if (is.null(x)) next
        pp[,count] <- p[,pt] + x*(p[,npt] - p[,pt])
        count <- count + 1
       
      }
    }
    p <- pp[,1:(count-1)]
  }

  if (any(zs)){
    z <- p[dnum,]
    n <- length(z)
    idx <- c(1:n,1)
    ii <- intersectHeight(z,hmin)
    ip <- z-hmin >= 0
    pp <- matrix(NA, nr=3, nc=2*n)
    count <- 1
    pt <- 1
    for (pt in 1:n){
      if (ip[pt]){
        pp[,count] <- p[,pt]
        count <- count + 1
      }
      if (ii[pt]==1){
        npt <- idx[pt+1]
        x <- calcIntersect(z[pt], z[npt], hmin)
        if (is.null(x)) next
        pp[,count] <- p[,pt] + x*(p[,npt] - p[,pt])
        count <- count + 1
       
      }
    }
    p <- pp[,1:(count-1)]
  }
  
  return(pp[,1:(count-1)])
}
  

  

  
splitHeight <- function(mshf, h, dnum=2){

  nstrips <- length(h)-1


  strips <- list()
  node.names <- list()
  face.names <- list()
  areas <- list()
  normals <- list()
  
  for (i in 1:nstrips){
    hmin <- h[i]
    hmax <- h[i+1]
    faces <- list()
    nn <- list()
    aa <- list()
    nodes <- character(0)
    fnames <- character(0)
    for (fn in names(mshf)){
      f <- mshf[[fn]]
      for (nm in names(f)){
        v <- f[[nm]]
        count <- 1
        np <- length(v$poly)
        tmp <- list()
        nntmp <- list()
        aatmp <- list()
        for (p  in v$poly){
          pp <- projectIntoHeight(p, hmin, hmax, dnum)
          if (!is.null(pp)){
            A <- poly3dArea(pp)
            if (A>=0){
              tmp[[count]] <- pp
              nntmp[[count]] <- poly3dNorm(pp)
              aatmp[[count]] <- A
              count <- count + 1
            }
          }
        }
        if (count > 1) nodes <- c(nodes, rep(nm, count-1))
        faces <- c(faces, tmp)
        fnames <- c(fnames, rep(fn, length(tmp)))
        aa <- c(aa, aatmp)
        nn <- c(nn, nntmp)
        
      }
      
      strips[[i]] <- faces
      node.names[[i]] <- nodes
      face.names[[i]] <- fnames
      areas[[i]] <- sapply(aa, function(x) x)
      normals[[i]] <- sapply(nn, function(x) x)
    }
  }
  return(list(faces=strips, nodes=node.names, fnames=face.names, h=h, areas=areas, normals=normals))
}

  

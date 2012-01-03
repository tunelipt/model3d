
polygon2tri <- function(p){

  d <- dim(p)

  np <- d[2]
  nd <- d[1]

  ntri <- np-2
  
  if (np==3){
    tri <- array(p, c(nd, 3, 1), list(c('x', 'y', 'z'), c('V1', 'V2', 'V3'), NULL))
    class(tri) <- 'mesh'
    return(tri)
  }
  xyz <- c('x', 'y', 'z')
  tri <- array(p, c(nd, 3, ntri), list(xyz[1:nd], c('V1', 'V2', 'V3'), NULL))

  for (i in 1:ntri)
    tri[,,i] <- p[,c(1, i+1, i+2)]

  if (nd==2)
    class(tri) <- 'mesh'
  else
    class(tri) <- 'mesh3d'
  
  return(tri)
  
    
  }

  
pmesh <- function(plst){

  if ('polygon3d' %in% class(plst[[1]]))
    class(plst) <- 'pmesh3d'
  else
    class(plst) <- 'pmesh'
  
  return(pmesh)
}


meshArea <- function(m) UseMethod("meshArea")

meshArea.mesh <- function(m) apply(m, 3, polygonArea.polygon)
meshArea.mesh3d <- function(m) apply(m, 3, polygonArea.polygon3d)
meshArea.pmesh <- function(m) sapply(m, function(p) polygonArea.polygon(p))
meshArea.pmesh3d <- function(m) sapply(m, function(p) polygonArea.polygon3d(p))


meshNormal <- function(m) UseMethod('meshNormal')

meshNormal.mesh <- function(m) rbind(x=rep(0, dim(m)[3]), y=rep(0, dim(m)[3]), z=rep(1, dim(m)[3]))
meshNormal.pmesh <- function(m) rbind(x=rep(0, length(m)), y=rep(0, length(m)), z=rep(1, length(m)))

meshNormal.mesh3d <- function(m) apply(m, 3, polygonNormal.polygon3d)

meshNormal.pmesh3d <- function(m) sapply(m, polygonNormal.polygon3d)


mergeMeshList <- function(trilst){

  if (!is.list(trilst)) trilst <- list(trilst)

  n <- length(trilst)
  ntri <- 0
  for (i in 1:n) ntri <- ntri + dim(trilst[[i]])[3]

  x <- array(dim=c(3, 3, ntri))

  count <- 1
  for (i in 1:n)
    for (k in 1:dim(trilst[[i]])[3]){
      x[,,count] <- trilst[[i]][,,k]
      count <- count + 1
    }
  class(x) <- class(trilst[[1]])
  return(x)
}  

pmesh2mesh <- function(pm){

  tri <- lapply(pm, polygon2tri)

  return(mergeMeshList(tri))
}

mesh2pmesh <- function(m){

  d <- dim(m)
  if (d[1]==2){
    cm <- "pmesh"
    cp <- 'polygon'
  }else{
    cm <- "pmesh3d"
    cp <- 'polygon3d'
  }
  
  ntri <- d[3]

  pm <- lapply(1:ntri, function(i){p <- m[,,i]; class(p) <- cp; p})
  class(pm) <- cm

  return(pm)
}
                                   
  
    
meshMeanNormal <- function(m){
  nn <- meshNormal(m)
  aa <- meshArea(m)
  A <- rbind(aa, aa, aa)
  return(rowSums(nn*A) / sum(aa))
}




                
meshSize <- function(m) UseMethod("meshSize")

meshSize.mesh <- function(m) dim(m)[3]
meshSize.pmesh <- function(m) length(m)
meshSize.mesh3d <- function(m) dim(m)[3]
meshSize.pmesh3d <- function(m) length(m)


meshGet <- function(m,i) UseMethod("meshGet")
meshGet.mesh <- function(m,i){p <- m[,,i]; class(p) <- 'polygon'; return(p)}
meshGet.pmesh <- function(m,i) m[[i]]
meshGet.mesh3d <- function(m,i) {p <- m[,,i]; class(p) <- 'polygon3d'; return(p)}
meshGet.pmesh3d <- function(m,i) m[[i]]


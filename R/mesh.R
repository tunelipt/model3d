#' Decomposes a polygon in triangles.
#'
#' Given a a simple convex polygon, this function decomposes it into triangles. The result is a mesh object (mesh3d for 3d polygons) that is a 3D array where the first index is the coordinate (x, y or z), the second the vertex (1, 2, 3) and the third the triangle number.
#' 
#' @param p Polygon
#' @return Array containing the triangles.
#' @export
#' @examples
#' p <- newPolygon(c(0, 1, 2, 1, 0), c(0, 0, 1, 2, 2))
#' tri <- polygon2tri(p)
#' plot(tri[1,,], tri[2,,])
#' for (i in 1:dim(tri)[3]) polygon(tri[1,,i], tri[2,,i])
polygon2tri <- function(p) UseMethod("polygon2tri")


#' Decomposes a polygon in triangles.
#'
#' Method for polygon class.
#' @export
polygon2tri.polygon <- function(p){

  d <- dim(p)

  np <- d[2]

  ntri <- np-2
  
  xy <- c('x', 'y')
  tri <- array(p, c(2, 3, ntri), list(xy, c('V1', 'V2', 'V3'), NULL))

  for (i in 1:ntri)
    tri[,,i] <- p[,c(1, i+1, i+2)]

  class(tri) <- 'mesh'
  
  return(tri)
  
    
  }


#' Decomposes a polygon in triangles.
#'
#' Method for polygon3d class.
#' @export
polygon2tri.polygon3d <- function(p){

  d <- dim(p)

  np <- d[2]

  ntri <- np-2
  
  xyz <- c('x', 'y', 'z')
  tri <- array(p, c(3, 3, ntri), list(xyz, c('V1', 'V2', 'V3'), NULL))
  
  for (i in 1:ntri)
    tri[,,i] <- p[,c(1, i+1, i+2)]
  
  class(tri) <- 'mesh3d'
  
  return(tri)
}


#' Creates a mesh from a list of polygons.
#'
#' This function basically gives the list the approppriate class name (pmesh or pmesh3d) depending on the class of the first element of the list.
#' @export
pmesh <- function(plst){

  if ('polygon3d' %in% class(plst[[1]]))
    class(plst) <- 'pmesh3d'
  else
    class(plst) <- 'pmesh'
  
  return(pmesh)
}


#' Compute the area of a mesh.
#'
#' A generic function that computes the area of a mesh.
#' @export
meshArea <- function(m) UseMethod("meshArea")

#' Compute the area of a mesh.
#'
#' Method for mesh.
#' @export
meshArea.mesh <- function(m) apply(m, 3, polygonArea.polygon)

#' Compute the area of a mesh.
#'
#' Method for mesh3d.
#' @export
meshArea.mesh3d <- function(m) apply(m, 3, polygonArea.polygon3d)


#' Compute the area of a mesh.
#'
#' Method for pmesh.
#' @export
meshArea.pmesh <- function(m) sapply(m, function(p) polygonArea.polygon(p))

#' Compute the area of a mesh.
#'
#' Method for pmesh3d.
#' @export
meshArea.pmesh3d <- function(m) sapply(m, function(p) polygonArea.polygon3d(p))


#' Compute the normal of each polygon in a mesh.
#'
#' Compute the normal of each polygon in a mesh. The normals are returned as a mtrix where each column correspondes to a polygon.
#' @export
meshNormal <- function(m) UseMethod('meshNormal')

#' Compute the normal of each polygon in a mesh.
#'
#' Method for mesh.
#' @export
meshNormal.mesh <- function(m) rbind(x=rep(0, dim(m)[3]), y=rep(0, dim(m)[3]), z=rep(1, dim(m)[3]))

#' Compute the normal of each polygon in a mesh.
#'
#' Method for pmesh.
#' @export
meshNormal.pmesh <- function(m) rbind(x=rep(0, length(m)), y=rep(0, length(m)), z=rep(1, length(m)))

#' Compute the normal of each polygon in a mesh.
#'
#' Method for mesh3d.
#' @export
meshNormal.mesh3d <- function(m) apply(m, 3, polygonNormal.polygon3d)

#' Compute the normal of each polygon in a mesh.
#'
#' Method for pmesh3d.
#' @export
meshNormal.pmesh3d <- function(m) sapply(m, polygonNormal.polygon3d)


#' Merges a list of triangular meshes.
#'
#' This function merges a list of triangular meshes into a single mesh.
#' @export
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

#' Converts a pmesh to mesh.
#'
#' Given a list of polygons, this function converts the polygons to triangles and merges them together to form a single triangular mesh.
#' @export
pmesh2mesh <- function(pm){

  tri <- lapply(pm, polygon2tri)

  return(mergeMeshList(tri))
}

#' Converts a mesh into a pmesh and mesh3d into pmesh3d.
#' @export
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
                                   
  
#' Mean normal on a surface.
#'
#' Calculates the mean normal of a set of polygons, the mean is weighed by the area of the polygons.
#' @export
meshMeanNormal <- function(m){
  nn <- meshNormal(m)
  aa <- meshArea(m)
  A <- rbind(aa, aa, aa)
  nm <- rowSums(nn*A) / sum(aa)
  return(nm / sqrt(sum(nm*nm)))
}





#' Number of polygons in a mesh.
#' @export
meshSize <- function(m) UseMethod("meshSize")

#' Number of polygons in a mesh.
#' @export
meshSize.mesh <- function(m) dim(m)[3]

#' Number of polygons in a pmesh.
#' @export
meshSize.pmesh <- function(m) length(m)

#' Number of polygons in a mesh3d.
#' @export
meshSize.mesh3d <- function(m) dim(m)[3]

#' Number of polygons in a pmesh3d.
#' @export
meshSize.pmesh3d <- function(m) length(m)


#' Return the ith polygon of a mesh.
#' @export
meshGet <- function(m,i) UseMethod("meshGet")

#' Return the ith polygon of a mesh.
#' @export
meshGet.mesh <- function(m,i){p <- m[,,i]; class(p) <- 'polygon'; return(p)}

#' Return the ith polygon of a pmesh.
#' @export
meshGet.pmesh <- function(m,i) m[[i]]

#' Return the ith polygon of a mesh3d.
#' @export
meshGet.mesh3d <- function(m,i) {p <- m[,,i]; class(p) <- 'polygon3d'; return(p)}

#' Return the ith polygon of a pmesh3d.
#' @export
meshGet.pmesh3d <- function(m,i) m[[i]]



#' Concatenates a list of lists into a large list.
#' @export
joinPolygons <- function(plst){

  p <- list()
  for (pp in plst)
    p <- c(p, pp)
  if (length(plst) > 0){
    class(p) <- class(plst[[1]])
    return(p)
  }
  return(p)
}

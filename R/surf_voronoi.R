#' Projection information of a face.
#'
#' Given a plane of projection identified by its orthonormal basis,
#' this function uses \code{\link{projectionInfo}} to calculate the projection
#' information so that points lying on the projection can be projected back to
#' 3D.
#'
#' @param face Polygonal mesh.
#' @param Orthonormal basis of projection
#' @param facor Whether the direction of normals should be inverted.
#' @return List of projection info for each polygon in the mesh.
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

  
#' Computes the voronoi diagram of points in a curved surface.
#'
#' Computing voronoi diagrams for 2D points is simple if libraries are
#' available. If the points are contrained in a curved 3D surface, things
#' are a little more complicated. This function projects the surface on a plane
#' computes the voronoi diagram in this plane and then projects the resulting
#' polygons back to 3D.
#'
#' @param m Polygonal 3D mesh.
#' @param x X Coordinate of the points that will be used to compute the Voronoi diagram
#' @param y Y Coordinate of the points that will be used to compute the Voronoi diagram
#' @param z Z Coordinate of the points that will be used to compute the Voronoi diagram
#' @param pnames Names of each point.
#' @param basefun Function used to obtain the projection plane.
#' @param ptsnames Names of the points.
#' @return A list containing 3D polygonal meshes of the influence area for each node.
#' @export
surfaceVoronoi <- function(m,  pts, basefun=getBase, ptsnames=NULL){

  
  x <- pts$x
  y <- pts$y
  z <- pts$z

  
  mean.norm <- meshMeanNormal(m)

  base2d <- basefun(mean.norm)[,1:2]

  pinfo <- faceProjectionInfo(m, base2d)
  

  pts <- pointsProject(rbind(x,y,z), base2d)
  npts <- length(x)
  if (npts < 1) return(NULL)
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
      pi3d <- project2dTo3d(pi[[1]], pinfo[[k]])
      tmp[[count]] <- pi3d
    }
    class(tmp) <- 'pmesh3d'
    plst[[i]] <- tmp
  }
  if (!is.null(ptsnames)) names(plst) <- ptsnames
  
  return(plst)
}

 

    
    

  

  
  
  

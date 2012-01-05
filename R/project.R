#' Get a 3D basis from a normal
#'
#' When projecting a 3D object into a plane it is necessary to have a
#' local 2D coordinate system on the plane. This coordinate system is given
#' by 2 orthonormal unit vectors. On the other hand, the plane o projection
#' is usually characterized by its normal so it is useful to get the local
#' 2D coordinate system from this normal. There are infinite solutions to
#' this problem but what this function do is use the normal (which in the
#' local reference system is z) to build "nice" local coordinate system.
#' This function does that by trying to align the local system as closely
#' as possible to (x, y, z) global coordinate system.
#'
#' @param bz A normal to the plane of projection.
#' @return Matrix where the columns are the unit vectors of the local coordinate system (x, y and z)
#' @export
#' @examples
#' bz <- c(sqrt(2)/2, sqrt(2,2), 0)
#' print(getBase(bz))
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

#' Projects 3D points into a plane.
#'
#' 3D points arranged as a matrix where each column corresponds
#' to a point is projected into a plane that has to direction unit vectors
#' x and y given by the columns of \code{base2d}.
#'
#' @param p Points
#' @param base2d Local plane coordinate system
#' @return a matrix where each column is a projected 2D point.
#' @export
pointsProject <- function(p, base2d) t(base2d) %*% p




#' Project a 3D polygon into a plane.
#'
#' Projects all vertices of a 3D polygon into a plane characterized by
#' 2 unit direction vectors.
#'
#' @param p 3D polygon to be projected.
#' @param base2d Matrix where each column corresponds to the local x and y dirs.
#' @return A 2D polygon.
#' @export
polygonProject <- function(p, base2d){
  p2 <- pointsProject(p, base2d)
  class(p2) <- 'polygon'
  return(p2)
}



#' Project a 3D face into a plane.
#'
#' Generic function that projects a 3D face (a polygon or a mesh) into a
#' plane characterized by 2 orhtonormal vectors.
#'
#' @param f Face
#' @param base2d Orthonormal local basis
#' @return Projection of f into base2d
#' @export
faceProject2d <- function(f, base2d) UseMethod('faceProject2d')

#' Method that projects a 3D points  into a plane.
#'
#' @param f Matrix where each row is a 3D point
#' @param base2d Orthonormal local basis
#' @return 2D projection
#' @export
faceProject2d.default <- function(f, base2d) pointsProject(f, base2d)

#' Method that projects a 3D polygon  into a plane.
#'
#' @param f polygon3d instance
#' @param base2d Orthonormal local basis
#' @return 2D projection, polygon
#' @export
faceProject2d.polygon3d <- function(f, base2d) polygonProject(f, base2d)


#' Method that projects a 3D triangular mesh  into a plane.
#'
#' @param f mesh3d instance
#' @param base2d Orthonormal local basis
#' @return 2D projection, mesh instance
#' @export
faceProject2d.mesh3d <- function(f, base2d){
  f2 <- apply(f, c(2,3), function(x) pointsProject(x, base2d))
  class(f2) <- 'mesh'
  return(f2)
}

#' Method that projects a 3D polygon mesh  into a plane.
#'
#' @param f pmesh3d instance
#' @param base2d Orthonormal local basis
#' @return 2D projection, pmesh instance
#' @export
faceProject2d.pmesh3d <- function(f, base2d){
  f2 <- lapply(f, function(p) faceProject2d.polygon3d(p, base2d))
  class(f2) <- 'pmesh'
  return(f2)
}


#' Projection information used to retrieve 3D data from 2D projections.
#'
#' When points that belong in a plane are projected into another plane,
#' some information is lost. This function gathers enough information from
#' the plane thar is being projected so that given any point *belonging* to
#' this plane can be determined from its 2D projection.
#'
#' @param p A point belonging to the plane that is being projected
#' @param nn A normal to the plane that is being projected
#' @param based2d Orthonormal basis of the projection plane.
#' @return Data structure that can be used to retrieve 3D data from its 2D projection.
#' @export
projectionInfo <- function(p, nn, base2d){

  A <- solve(rbind(t(base2d), nn))  
  k <- sum(p*nn)
  pinfo <- list(p=p, norm=nn, A=A, base2d=base2d, k=k)
  class(pinfo) <- 'projectionInfo'
  return(pinfo)
  
}


#' Project a 3D face into a plane.
#'
#' Generic function that projects a 3D face (a polygon or a mesh) into a
#' plane characterized by 2 orhtonormal vectors.
#'
#' @param f Face
#' @param base2d Orthonormal local basis
#' @return Projection of f into base2d
#' @seealso \code{\link{faceProject2d}}
#' @export
project3dTo2d <- function(p, base2d) UseMethod('project3dTo2d')

#' Method that projects a 3D points  into a plane.
#'
#' @param f Matrix where each row is a 3D point
#' @param base2d Orthonormal local basis
#' @return 2D projection
#' @seealso \code{\link{faceProject2d}}, \code{\link{pointsProject}}
#' @export
project3dTo2d.default <- function(p, base2d) faceProject2d.default(p, base2d)


#' Method that projects a 3D polygon  into a plane.
#'
#' @param f polygon3d instance
#' @param base2d Orthonormal local basis
#' @return 2D projection, polygon
#' @seealso \code{\link{faceProject2d}}, \code{\link{polygonProject}}
#' @export
project3dTo2d.polygon3d <- function(p, base2d) faceProject2d.polygon3d(p, base2d)


#' Method that projects a 3D triangular mesh  into a plane.
#'
#' @param f mesh3d instance
#' @param base2d Orthonormal local basis
#' @return 2D projection, mesh instance
#' @export
project3dTo2d.mesh3d <- function(p, base2d) faceProject2d.mesh3d(p, base2d)

#' Method that projects a 3D polygon mesh  into a plane.
#'
#' @param f pmesh3d instance
#' @param base2d Orthonormal local basis
#' @return 2D projection, pmesh instance
#' @export
project3dTo2d.pmesh3d <- function(p, base2d) faceProject2d.pmesh3d(p, base2d)


#' Get back 3D data from 2D projection.
#'
#' Retrieves initial 3D data from the 2D projection.
#'
#' @param p 2D projection
#' @param pinfo Projection info obtained from \code{\link{projectionInfo}}
#' @return 3D information.
#' @seealso \code{\link{projectionInfo}}, \code{\link{project3dTo2d}}
#' @export
project2dTo3d <- function(p, pinfo) UseMethod('project2dTo3d')

#' Get back 3D points from 2D projection.
#'
#' Retrieves initial 3D points from ist 2D projection.
#'
#' @param p 2D projection of points
#' @param pinfo Projection info obtained from \code{\link{projectionInfo}}
#' @return 3D points.
#' @seealso \code{\link{projectionInfo}} , \code{\link{project3dTo2d}}
#' @export
project2dTo3d.default <- function(p, pinfo) pinfo$A %*% c(p, pinfo$k)


#' Get back 3D polygon from 2D projection.
#'
#' Retrieves initial 3D polygon from ist 2D projection.
#'
#' @param p 2D projection of polygon
#' @param pinfo Projection info obtained from \code{\link{projectionInfo}}
#' @return 3D polygon.
#' @seealso \code{\link{projectionInfo}} , \code{\link{project3dTo2d}}
#' @export
project2dTo3d.polygon <- function(p, pinfo){
  p3 <- pinfo$A %*% rbind(p, pinfo$k)
  class(p3) <- 'polygon3d'
  return(p3)
}

#' Get back 3D polygon from 2D projection.
#'
#' Retrieves initial 3D polygon from ist 2D projection.
#'
#' @param p 2D projection of polygon, gpc.poly format.
#' @param pinfo Projection info obtained from \code{\link{projectionInfo}}
#' @return 3D polygon.
#' @seealso \code{\link{projectionInfo}} , \code{\link{project3dTo2d}}
#' @export
project2dTo3d.gpc.poly <- function(p, pinfo) project2dTo3d.polygon(gpc2poly(p), pinfo)



  

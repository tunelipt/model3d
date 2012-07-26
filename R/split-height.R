

#' Calculates where the segments of a polygon intersect a certain height.
#'
#' Given a polygone, if one wants to chop it at a given height, it is
#' important to figure out which segments of the polygon cross this height.
#' This function identifies which segments cross the height. One point
#' should be remembered: sometimes a vertex has height h. In this case
#' it should be handled differently. Edges that cross the height h are
#' identified by flag 1. If the first vertex of the edge intersects the height,
#' the edge is identified by flag 2. Otherwise the edge is identifed by flag 0
#'
#' @param z Heights of each vertex of the polygon.
#' @param h height with which to cut the polygon.
#' @return Integer vector containing intersection information.
#' @export
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

#' Calculates the position of a point along a segment.
#'
#' Thus function calculates the position along a segment where a point
#' is located. The position is non-dimensional in the sense that the
#' first point is at position 0 and the second point is at position 1.
#'
#' @param z1 First point.
#' @param z2 Second point.
#' @param z Point whose position is to be determined.
#' @return The linear position of the point where z1 is 0 and z2 is 1.
#' @export
calcIntersect <- function(z1, z2, z){

  dz <- z2-z1

  if (dz==0) return(NULL)

  return( (z-z1)/dz )
}


  
#' Chops a polygon at two different heights.
#'
#' Given a polygon and two levels, this function chops the polygon,
#' returning a polygon that is the intersection between the original polygon
#' and a strip from the lower height to the upper height. If there is not
#' intersection this function returns \code{NULL}.
#' The height can be any coordinate axis: x -> 1, y -> 2 and z -> 3 (default).
#'
#' @param p Polygon, a matrix where the rows are x, y, z coordinates of the vertices.
#' @param hmin Minimum height to chop the polygon.
#' @param hmax Maximum height to chop the polygon.
#' @param dnum Axis along which the heights are defined.
#' @return Chopped polygon or NULL if there is no intersection.
#' @export
projectIntoHeight <- function(p, hmin, hmax, dnum=3){
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
  

#' Splits a meshed surface in strips for a set of heights.
#'
#' This function splits a face along strips defined by
#' positions along a coordinate axis.
#'
#' @param List of node influence area meshes.
#' @param h Heights to be used to split the face.
#' @param If more than one face should be splitted this argument accumulates the result.
#' @param dnum Direction used to split x->1, y->2 and z->3.
#' @return A list of node influence area for each height strip.
#' @export
splitHeight <- function(face, h, mshh=NULL, dnum=3){

  nstrips <- length(h)-1
  
  if (is.null(mshh)){
    mshh <- list()
    for (i in 1:nstrips)
      mshh[[i]] <- list()
    attributes(mshh) <- list(h=h)
  }


  for (i in 1:nstrips){
    hmin <- h[i]
    hmax <- h[i+1]
    for (nname in names(face)){
      node <- face[[nname]]
      for (p in node){
        hlim <- range(p[dnum,])
        if (hmax <= hlim[1] || hmin >= hlim[2]) next
        pchop <- projectIntoHeight(p, hmin, hmax, dnum)
        if (is.null(mshh[[i]][[nname]]))
          mshh[[i]][[nname]] <- list(pchop)
        else
          mshh[[i]][[nname]] <- append(mshh[[i]][[nname]], list(pchop))
      }
    }
  }

  for (i in 1:nstrips)
    if (length(mshh[[i]]) > 0)
      for (k in 1:length(mshh[[i]]))
        class(mshh[[i]][[k]]) <- 'pmesh3d'
      
  
  return(mshh)
      
}

#' Splits all faces of a building along a direction.
#'
#' This function splits all faces of a building or structure along a coordinate axis.
#' It simply calls, repeatedly the function \code{\link{splitHeight}} for multiple faces.
#'
#' @param flst List of discretizations of each face.
#' @param h Positions defining where the faces should be split.
#' @param dnum Direction where the faces should be split. x->1, y->2, z->3.
#' @return Meshes of each influence area for each level.
#' @export
splitBuilding <- function(flst, h, dnum=3){

  mshh <- NULL

  for (f in flst)
    mshh <- splitHeight(f, h, mshh, dnum)

  return(mshh)
}



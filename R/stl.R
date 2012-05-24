



#' Reads an STL file
#'
#' Reads a binary STL format into a mesh.
#'
#' @param fname Name of stl file.
#' @return A \code{mesh3d} mesh.
#' @export
stlRead <- function(fname){


  con <- file(fname, open='rb')
  on.exit(close(con))
  header <- readBin(con, what=raw(1), n=80, size=1)

  nfacets <- readBin(con, what=integer(1), n=1, size=4)

  normals <- matrix(0.0, nr=3, nc=nfacets, dimnames=list(c('x', 'y', 'z'), NULL))
  mesh <- array(0.0, dim=c(3,3,nfacets), dimnames=list(c('x', 'y', 'z'), c('V1', 'V2', 'V3'), NULL))

  for (i in 1:nfacets){
    ff <- readBin(con, what=double(1), n=12, size=4)
    abc <- readBin(con, what=integer(1), n=1, size=2, signed=FALSE)

    mesh[,,i] <- ff[4:12]
    normals[,i] <- ff[1:3]
  }

  class(mesh) <- 'mesh3d'
  return(mesh)
  #return(list(mesh=mesh, normal=normals))
}
    

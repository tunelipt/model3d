\name{faceProjectionInfo}
\alias{faceProjectionInfo}
\title{Projection information of a face.}
\usage{
  faceProjectionInfo(face, base2d, facor = 1)
}
\arguments{
  \item{face}{Polygonal mesh.}

  \item{Orthonormal}{basis of projection}

  \item{facor}{Whether the direction of normals should be
  inverted.}
}
\value{
  List of projection info for each polygon in the mesh.
}
\description{
  Given a plane of projection identified by its orthonormal
  basis, this function uses \code{\link{projectionInfo}} to
  calculate the projection information so that points lying
  on the projection can be projected back to 3D.
}


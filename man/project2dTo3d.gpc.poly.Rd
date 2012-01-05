\name{project2dTo3d.gpc.poly}
\alias{project2dTo3d.gpc.poly}
\title{Get back 3D polygon from 2D projection.}
\usage{
  project2dTo3d.gpc.poly(p, pinfo)
}
\arguments{
  \item{p}{2D projection of polygon, gpc.poly format.}

  \item{pinfo}{Projection info obtained from
  \code{\link{projectionInfo}}}
}
\value{
  3D polygon.
}
\description{
  Retrieves initial 3D polygon from ist 2D projection.
}
\seealso{
  \code{\link{projectionInfo}} ,
  \code{\link{project3dTo2d}}
}


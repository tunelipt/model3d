\name{splitBuilding}
\alias{splitBuilding}
\title{Splits all faces of a building along a direction.}
\usage{
  splitBuilding(flst, h, dnum = 3)
}
\arguments{
  \item{flst}{List of discretizations of each face.}

  \item{h}{Positions defining where the faces should be
  split.}

  \item{dnum}{Direction where the faces should be split.
  x->1, y->2, z->3.}
}
\value{
  Meshes of each influence area for each level.
}
\description{
  This function splits all faces of a building or structure
  along a coordinate axis. It simply calls, repeatedly the
  function \code{\link{splitHeight}} for multiple faces.
}


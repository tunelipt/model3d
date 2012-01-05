\name{getBase}
\alias{getBase}
\title{Get a 3D basis from a normal}
\usage{
  getBase(bz)
}
\arguments{
  \item{bz}{A normal to the plane of projection.}
}
\value{
  Matrix where the columns are the unit vectors of the
  local coordinate system (x, y and z)
}
\description{
  When projecting a 3D object into a plane it is necessary
  to have a local 2D coordinate system on the plane. This
  coordinate system is given by 2 orthonormal unit vectors.
  On the other hand, the plane o projection is usually
  characterized by its normal so it is useful to get the
  local 2D coordinate system from this normal. There are
  infinite solutions to this problem but what this function
  do is use the normal (which in the local reference system
  is z) to build "nice" local coordinate system. This
  function does that by trying to align the local system as
  closely as possible to (x, y, z) global coordinate
  system.
}
\examples{
bz <- c(sqrt(2)/2, sqrt(2,2), 0)
print(getBase(bz))
}


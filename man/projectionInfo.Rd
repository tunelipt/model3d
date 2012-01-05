\name{projectionInfo}
\alias{projectionInfo}
\title{Projection information used to retrieve 3D data from 2D projections.}
\usage{
  projectionInfo(p, nn, base2d)
}
\arguments{
  \item{p}{A point belonging to the plane that is being
  projected}

  \item{nn}{A normal to the plane that is being projected}

  \item{based2d}{Orthonormal basis of the projection
  plane.}
}
\value{
  Data structure that can be used to retrieve 3D data from
  its 2D projection.
}
\description{
  When points that belong in a plane are projected into
  another plane, some information is lost. This function
  gathers enough information from the plane thar is being
  projected so that given any point *belonging* to this
  plane can be determined from its 2D projection.
}


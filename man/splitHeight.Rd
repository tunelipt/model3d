\name{splitHeight}
\alias{splitHeight}
\title{Splits a meshed surface in strips for a set of heights.}
\usage{
  splitHeight(face, h, mshh = NULL, dnum = 3)
}
\arguments{
  \item{List}{of node influence area meshes.}

  \item{h}{Heights to be used to split the face.}

  \item{If}{more than one face should be splitted this
  argument accumulates the result.}

  \item{dnum}{Direction used to split x->1, y->2 and z->3.}
}
\value{
  A list of node influence area for each height strip.
}
\description{
  This function splits a face along strips defined by
  positions along a coordinate axis.
}


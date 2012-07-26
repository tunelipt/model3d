\name{intersectHeight}
\alias{intersectHeight}
\title{Calculates where the segments of a polygon intersect a certain height.}
\usage{
  intersectHeight(z, h)
}
\arguments{
  \item{z}{Heights of each vertex of the polygon.}

  \item{h}{height with which to cut the polygon.}
}
\value{
  Integer vector containing intersection information.
}
\description{
  Given a polygone, if one wants to chop it at a given
  height, it is important to figure out which segments of
  the polygon cross this height. This function identifies
  which segments cross the height. One point should be
  remembered: sometimes a vertex has height h. In this case
  it should be handled differently. Edges that cross the
  height h are identified by flag 1. If the first vertex of
  the edge intersects the height, the edge is identified by
  flag 2. Otherwise the edge is identifed by flag 0
}


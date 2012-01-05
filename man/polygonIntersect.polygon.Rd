\name{polygonIntersect.polygon}
\alias{polygonIntersect.polygon}
\title{Polygon Intersection.}
\usage{
  polygonIntersect.polygon(p1, p2)
}
\description{
  Method that computes intersection of two polygons of
  class polygon. It converts the polygons to gpc.poly
  polygons and then uses gpclib to compute the
  intersection. The result is then converted to a list of
  polygons of class polygon.
}


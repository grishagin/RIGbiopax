% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_deadend_refs.R
\name{remove_deadend_refs}
\alias{remove_deadend_refs}
\title{Remove Dead-End References from BioPAX}
\usage{
remove_deadend_refs(biopax)
}
\arguments{
\item{biopax}{A BioPAX object.}
}
\description{
Takes in a BioPAX file, and removes those instances that refer 
to instances that do not exist.
}
\author{
Ivan Grishagin
}

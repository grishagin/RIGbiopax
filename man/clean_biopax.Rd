% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_biopax.R
\name{clean_biopax}
\alias{clean_biopax}
\title{Clean BioPAX Object}
\usage{
clean_biopax(biopax)
}
\arguments{
\item{biopax}{BioPAX object.}
}
\description{
1) Strips hashes from property_value column.
2) Cleans a BioPAX object from html tags and garbled symbols. 
3) Unifies representation of some classes, common properties and property values.
4) Removes deadend references.
}
\details{

}
\author{
Ivan Grishagin
}

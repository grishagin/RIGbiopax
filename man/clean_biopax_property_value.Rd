\name{clean_biopax_property_value}
\alias{clean_biopax_property_value}
\title{
Clean up BioPAX property_value Column
}
\description{
Function removes some undesirable elements from BioPAX property_value column (see Details)
}
\usage{
clean_biopax_property_value(biopax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{biopax}{
  BioPAX object to be cleaned.
  }
}
\details{
-- removes html tags\cr 
-- removes eol chars\cr
-- replaces entire entries that start with "Recname:" with the recommended names
}
\value{
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ivan Grishagin
}
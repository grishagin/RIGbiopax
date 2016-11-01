\name{splitbiopax_from_dt}
\alias{splitbiopax_from_dt}
\title{
Convert a BioPAX-Style Data Table to List of BioPAX Objects (by Pathway)
}
\description{
Converts BioPAX-style data table to a list of BioPAX objects by pathway. Has an option of recording all of these BioPAX objects to file.
}
\usage{
splitbiopax_from_dt(dTable
                    ,write_to_files=FALSE
                    ,pwtoextract=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dTable}{
  BioPAX-style data table or data frame. 
  }
  \item{write_to_files}{
  logical. Write to files? Default is \code{FALSE}.
  }
  \item{pwtoextract}{
  vector of pathway ids. If \code{NULL} (default), all pathways are extracted.
  }
}
\details{
Returns a list of biopax objects.
}
\value{

}
\references{
}
\author{
Ivan Grishagin
}
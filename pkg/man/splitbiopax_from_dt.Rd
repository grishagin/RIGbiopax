\name{splitbiopax_from_dt}
\alias{splitbiopax_from_dt}
\title{
Convert BioPAX Object to List of BioPAX Objects (by Pathway)
}
\description{
Splits a BioPAX object into a list of BioPAX objects by pathway. Has an option of recording all of these BioPAX objects to files.
}
\usage{
splitbiopax_from_dt(dTable
                    ,write_to_files=FALSE
                    ,pwtoextract=NULL
                    ,unifyids=TRUE
                    )
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
  \item{unifyids}{
  logical. Unify IDs by making all elements of one class receive ids in the form "ClassNUMBER", with NUMBER ranging from one to number of elements of such class?
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
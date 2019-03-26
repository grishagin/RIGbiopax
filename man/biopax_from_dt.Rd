\name{biopax_from_dt}
\alias{biopax_from_dt}
\title{
Convert a BioPAX-Style Data Table to BioPAX Object
}
\description{
Converts BioPAX-style data table to a BioPAX object. Has an option of recording the BioPAX object to file.
}
\usage{
biopax_from_dt(dTable
               ,filename=NULL
               ,encoding="ISO-8859-1")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dTable}{
  BioPAX-style data table.
  }
  \item{filename}{
  output filename. If \code{NULL} (default), does not output to file.
  }
    \item{encoding}{
  output XML file encoding.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ivan Grishagin
}
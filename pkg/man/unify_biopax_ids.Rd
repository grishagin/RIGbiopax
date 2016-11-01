\name{unify_biopax_ids}
\alias{unify_biopax_ids}
\title{
Convert a BioPAX-Style Data Table to BioPAX Object
}
\description{
Converts ids in a BioPAX-style data table to a uniform format \code{CLASS#}, where \code{CLASS} is an element class, and \code{#} is a number.
}
\usage{
unify_biopax_ids(biopax_dt
                 ,idtag=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{biopax_dt}{
  BioPAX-style data table.
  }
\item{idtag}{
  string to append to each id (defaults to \code{NULL}).
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
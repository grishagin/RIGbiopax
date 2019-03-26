\name{writeBiopax_Rancho}
\alias{writeBiopax_Rancho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Write BioPAX Object to File
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Writes a BioPAX object to the XML (\code{*.owl}) file. Standard \code{writeBiopax} function from package \code{rBiopaxParser} has poor performance that is felt for anything larger than 10,000 lines.
}
\usage{
writeBiopax_Rancho(biopax
                  ,filename = NULL
                  ,overwrite = FALSE
                  ,biopaxlevel=3
                  ,namespaces=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{biopax}{
  biopax object.
  }
  \item{filename}{
  desired filename to produce. If \code{NULL}, the name will be generated automatically.
  }
    \item{overwrite}{
  logical. Overwrite?
  }
    \item{biopaxlevel}{
  either 2 or 3. Affects the namespaces.
  }
    \item{namespaces}{
  just in case custom namespaces are desired.
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
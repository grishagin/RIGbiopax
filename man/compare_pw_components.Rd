\name{compare_pw_components}
\alias{compare_pw_components}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compare Proteins Extracted from BioPAX Pathway to Source Table
}
\description{
Internal function. For a given pathway, perform comparison of proteins extracted from BioPAX to those listed in the source table/dataframe.
}
\usage{
compare_pw_components(biopax_prot = NULL
                      ,toxdb_prot = NULL
                      ,dfrow = NULL
                      ,output = c("all", "match", "biopax", "toxdb"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{biopax_prot}{
%%     ~~Describe \code{biopax_prot} here~~
}
  \item{toxdb_prot}{
%%     ~~Describe \code{toxdb_prot} here~~
}
  \item{dfrow}{
%%     ~~Describe \code{dfrow} here~~
}
  \item{output}{
%%     ~~Describe \code{output} here~~
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
\author{
Ivan Grishagin
}
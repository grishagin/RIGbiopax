\name{expand.df.via.split.col}
\alias{expand.df.via.split.col}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Split Column(s) Based on Pattern
}
\description{
Wrapper for a \code{strsplit} function. In a given dataframe, split specified column values based on certain pattern. New values are inserted by duplicating corresponding rows, where the cells that were affected were located. Can work for multiple columns.
}
\usage{
expand.df.via.split.col(dFrame = NULL
                        ,colsToSplit = NULL
                        ,patternToSplit = ",|;|\\|")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dFrame}{
%%     ~~Describe \code{dFrame} here~~
}
  \item{colsToSplit}{
%%     ~~Describe \code{colToSplit} here~~
}
  \item{patternToSplit}{
%%     ~~Describe \code{patternToSplit} here~~
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
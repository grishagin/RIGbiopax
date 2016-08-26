\name{expand.df.via.split.col}
\alias{expand.df.via.split.col}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Split Column(s) Based on Pattern
}
\description{
Wrapper for a \code{strsplit} function. In a given dataframe, split specified column values based on certain pattern. 
With multiple columns, it can work in two ways: either one after the other with \code{at_once=FALSE} or all at the same time with \code{at_once=TRUE} (default). 
When all columns are processed at once, the function attempts to find the same number of elements (that are produced after splitting) in all columns for every row that contains the pattern.
 If some columns have just one element, it is replicated to match the maximum number of elements. If the number of elements is otherwise not equal between the columns, an error is thrown.
}
\usage{
expand.df.via.split.col(dFrame
                        ,colsToSplit
                        ,patternToSplit = "\\\\|"
                        ,at_once=TRUE)
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
  \item{at_once}{
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
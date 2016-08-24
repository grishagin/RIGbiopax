\name{shrink.df.via.merge.col}
\alias{shrink.df.via.merge.col}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Merge Columns in Dataframe Using Pattern
}
\description{
Looks up duplicate values in the \code{colKey} column, and then merges all corresponding values in the \code{colToMerge} column, separating them with a \code{patternToMerge} separator.
}
\usage{
shrink.df.via.merge.col(dFrame = NULL
                        ,colKey = NULL
                        ,colToMerge = NULL
                        ,patternToMerge = ",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dFrame}{
%%     ~~Describe \code{dFrame} here~~
}
  \item{colKey}{
  key column. Duplicate values are looked up in this column.
}
  \item{colToMerge}{
  columns to merge and return. If \code{NULL}, merges returns all but \code{colKey} columns.
}
  \item{patternToMerge}{
  pattern to use for merging of the values. Defaults to comma (\code{,}).
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
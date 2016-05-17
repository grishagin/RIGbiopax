\name{expand.df.via.split.col}
\alias{expand.df.via.split.col}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Split Column Based on Pattern
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
expand.df.via.split.col(dFrame = NULL, colToSplit, patternToSplit = ",|;|\\|")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dFrame}{
%%     ~~Describe \code{dFrame} here~~
}
  \item{colToSplit}{
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dFrame = NULL, colToSplit, patternToSplit = ",|;|\\|") 
{
    rowsAffected <- grep(patternToSplit, dFrame[, colToSplit[1]])
    if (length(rowsAffected) > 0) {
        splitList <- dFrame[rowsAffected, colToSplit[1]] \%>\% 
            strsplit(split = patternToSplit)
        rowReplicates <- rep(1, nrow(dFrame))
        rowReplicates[rowsAffected] <- sapply(splitList, length)
        finalRowsToTake <- rep(1:nrow(dFrame), rowReplicates)
        dFrame <- dFrame[finalRowsToTake, ]
        rowsAffected_2 <- grep(patternToSplit, dFrame[, colToSplit[1]])
        dFrame[rowsAffected_2, colToSplit] <- unlist(splitList)
    }
    return(dFrame)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

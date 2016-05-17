\name{shrink.df.via.merge.col}
\alias{shrink.df.via.merge.col}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Merge Columns in Dataframe Using Pattern
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
shrink.df.via.merge.col(dFrame = NULL, colKey = NULL, colToMerge = NULL, patternToMerge = ",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dFrame}{
%%     ~~Describe \code{dFrame} here~~
}
  \item{colKey}{
%%     ~~Describe \code{colKey} here~~
}
  \item{colToMerge}{
%%     ~~Describe \code{colToMerge} here~~
}
  \item{patternToMerge}{
%%     ~~Describe \code{patternToMerge} here~~
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
function (dFrame = NULL, colKey = NULL, colToMerge = NULL, patternToMerge = ",") 
{
    duplicate_keys <- table(dFrame[, colKey]) \%>\% .[. > 1] \%>\% 
        names
    merged_id_df <- adply(.data = duplicate_keys, .margins = 1, 
        .fun = function(KEY) {
            merged_id_df <- dFrame[dFrame[, colKey] \%chin\% KEY, 
                colToMerge] \%>\% lapply(paste, collapse = patternToMerge) \%>\% 
                as.data.frame
            merged_id_df$KEY <- KEY
            merged_id_df <- merged_id_df \%>\% .[, c("KEY", colToMerge)]
            colnames(merged_id_df) <- c(colKey, colToMerge)
            return(merged_id_df)
        }) \%>\% dplyr::select(-1)
    dFrame <- dFrame \%>\% .[!(dFrame[, colKey] \%chin\% duplicate_keys), 
        c(colKey, colToMerge)] \%>\% filter(!(colKey \%chin\% duplicate_keys)) \%>\% 
        rbind(merged_id_df)
    return(dFrame)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

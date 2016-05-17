\name{add.biopax.ids}
\alias{add.biopax.ids}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
add.biopax.ids(all_pathways, pw_biopax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{all_pathways}{
%%     ~~Describe \code{all_pathways} here~~
}
  \item{pw_biopax}{
%%     ~~Describe \code{pw_biopax} here~~
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
function (all_pathways, pw_biopax) 
{
    all_pathways_pwid <- adply(.data = as.data.frame(all_pathways), 
        .margins = 1, .fun = function(dfrow) {
            biopaxID <- pw_biopax$biopax.Pathway.ID[pw_biopax$biopax.Pathway.Name \%chin\% 
                dfrow$biopax.Pathway.Name]
            dfrow <- dfrow \%>\% .[rep(1, length(biopaxID)), ] \%>\% 
                mutate(biopax.Pathway.ID = biopaxID)
            if (length(biopaxID) > 1) {
                message("Found ", length(biopaxID), " IDs for ", 
                  unique(dfrow$biopax.Pathway.Name))
            }
            return(dfrow)
        })
    return(all_pathways_pwid)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

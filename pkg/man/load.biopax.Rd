\name{load.biopax}
\alias{load.biopax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
load.biopax(source_name = NULL, source_dir = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{source_name}{
%%     ~~Describe \code{source_name} here~~
}
  \item{source_dir}{
%%     ~~Describe \code{source_dir} here~~
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
function (source_name = NULL, source_dir = NULL) 
{
    objects <- ls(envir = parent.env(env = environment()))
    is_present_biopax <- sapply(objects, function(obj) {
        "biopax" \%in\% class(get(obj))
    })
    if (sum(is_present_biopax) < 1) {
        biopax <- load.biopax.file(source_name = source_name, 
            source_dir = source_dir)
    }
    else if (sum(is_present_biopax) > 1) {
        stop("More than one biopax object present in environment already")
    }
    else {
        biopax <- get(subset(objects, is_present_biopax))
    }
    return(biopax)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

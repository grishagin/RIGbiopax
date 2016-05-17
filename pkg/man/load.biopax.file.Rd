\name{load.biopax.file}
\alias{load.biopax.file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Load BioPAX OWL File
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
load.biopax.file(source_name = NULL, source_dir = NULL)
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
    require(rBiopaxParser)
    if (is.null(source_name)) {
        stop("Provide source_name!")
    }
    if (is.null(source_dir)) {
        source_dir <- getwd()
    }
    owl.dirs <- list.dirs(source_dir, recursive = FALSE)
    owl_num <- agrep(source_name, owl.dirs, max.distance = 0.2)
    owl.filename <- list.files(owl.dirs[owl_num], pattern = "\\.owl", 
        full.names = TRUE)
    if (length(owl.filename) != 1) {
        stop("There's not just one biopax file in the following folder:\n", 
            owl.dirs[owl_num])
    }
    owl_biopax <- readBiopax(owl.filename)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{load.biopax.pathways}
\alias{load.biopax.pathways}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
load.biopax.pathways(owl_biopax = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{owl_biopax}{
%%     ~~Describe \code{owl_biopax} here~~
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
function (owl_biopax = NULL) 
{
    if (is.null(owl_biopax)) {
        stop("Biopax object was not provided")
    }
    name_properties <- grep("name", owl_biopax$dt$property, ignore.case = TRUE, 
        value = TRUE) \%>\% unique() \%>\% tolower
    if ("displayname" \%in\% name_properties) {
        owl_biopax_instances_stname <- selectInstances(biopax = owl_biopax, 
            class = "pathway", property = "displayName")
        pw_biopax <- owl_biopax_instances_stname \%>\% dplyr::select(id, 
            name = property_value)
    }
    else if ("standardname" \%in\% name_properties & "name" \%in\% 
        name_properties) {
        owl_biopax_instances_stname <- selectInstances(biopax = owl_biopax, 
            class = "pathway", property = "standardName")
        owl_biopax_instances_name <- selectInstances(biopax = owl_biopax, 
            class = "pathway", property = "name")
        matching_rows <- owl_biopax_instances_stname$id \%in\% 
            owl_biopax_instances_name$id
        pw_biopax <- owl_biopax_instances_stname[matching_rows, 
            ] \%>\% dplyr::select(id, name = property_value)
    }
    else {
        stop("load.biopax.pathways: check your name properties!")
    }
    pw_biopax$name <- tolower(pw_biopax$name)
    colnames(pw_biopax) <- c("biopax.Pathway.ID", "biopax.Pathway.Name")
    return(pw_biopax)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

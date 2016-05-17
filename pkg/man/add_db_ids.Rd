\name{add_db_ids}
\alias{add_db_ids}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
add_db_ids(owl_biopax = NULL, pw_id = NULL, biopaxlevel = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{owl_biopax}{
%%     ~~Describe \code{owl_biopax} here~~
}
  \item{pw_id}{
%%     ~~Describe \code{pw_id} here~~
}
  \item{biopaxlevel}{
%%     ~~Describe \code{biopaxlevel} here~~
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
function (owl_biopax = NULL, pw_id = NULL, biopaxlevel = 3) 
{
    if (is.null(owl_biopax) || is.null(pw_id)) {
        stop("Some of the required parameters were not provided.")
    }
    else if (is.na(pw_id)) {
        results_df <- data.frame(biopax.Component.ID = NA, biopax.Gene.ID.Type = NA, 
            biopax.Gene.ID = NA)
        return(results_df)
    }
    require(tidyr)
    require(dplyr)
    require(data.table)
    require(rBiopaxParser)
    pw_components_ids <- listPathwayComponents(biopax = owl_biopax, 
        id = pw_id, returnIDonly = TRUE, biopaxlevel = biopaxlevel)
    pw_components <- selectInstances(biopax = owl_biopax, id = pw_components_ids, 
        includeReferencedInstances = TRUE)
    print(pw_id)
    results_df <- pw_components \%>\% filter(class == "RelationshipXref" | 
        class == "UnificationXref") \%>\% dplyr::select(biopax.Component.ID = id, 
        property, property_value) \%>\% spread(property, property_value)
    if (nrow(results_df) < 1) {
        return(data.frame(biopax.Component.ID = NA, biopax.Gene.ID.Type = NA, 
            biopax.Gene.ID = NA))
    }
    setnames(results_df, "db", "biopax.Gene.ID.Type")
    setnames(results_df, "id", "biopax.Gene.ID")
    results_df$biopax.Gene.ID[results_df$biopax.Gene.ID \%chin\% 
        ""] <- NA
    results_df <- results_df \%>\% dplyr::select(biopax.Component.ID, 
        biopax.Gene.ID.Type, biopax.Gene.ID)
    return(results_df)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

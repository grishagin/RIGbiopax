\name{load.toxdb.genes.per.source}
\alias{load.toxdb.genes.per.source}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Load Dataframe with Genes/Pathways
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
load.toxdb.genes.per.source(source_name, source_dir = getwd(), toxdb_genes_file = "toxdb_pathways_15Jun.txt", all_pathways)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{source_name}{
%%     ~~Describe \code{source_name} here~~
}
  \item{source_dir}{
%%     ~~Describe \code{source_dir} here~~
}
  \item{toxdb_genes_file}{
%%     ~~Describe \code{toxdb_genes_file} here~~
}
  \item{all_pathways}{
%%     ~~Describe \code{all_pathways} here~~
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
function (source_name, source_dir = getwd(), toxdb_genes_file = "toxdb_pathways_15Jun.txt", 
    all_pathways) 
{
    options(stringsAsFactors = FALSE)
    toxdb <- read.table(file.path(source_dir, toxdb_genes_file), 
        header = TRUE, quote = "", sep = "\t") \%>\% dplyr::select(toxdb.Pathway.ID = Pathway.ID, 
        toxdb.Pathway.Name = Pathway.Name, toxdb.Gene.ID = Gene.ID, 
        toxdb.Gene.Symbol = Gene.Symbol)
    toxdb$toxdb.Pathway.Name <- tolower(toxdb$toxdb.Pathway.Name)
    toxdb <- toxdb \%>\% subset(.$toxdb.Pathway.ID \%in\% all_pathways$toxdb.Pathway.ID) \%>\% 
        mutate(Source = source_name)
    return(toxdb)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

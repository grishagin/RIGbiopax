\name{add.MULT.symbols.entrezids}
\alias{add.MULT.symbols.entrezids}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
add.MULT.symbols.entrezids(df_pw_proteins = NULL, keytypes = c("entrezgene", "mim", "uniprot", "unigene", "ensemblgene", "hprd", "hgnc"), filter_keytypes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df_pw_proteins}{
%%     ~~Describe \code{df_pw_proteins} here~~
}
  \item{keytypes}{
%%     ~~Describe \code{keytypes} here~~
}
  \item{filter_keytypes}{
%%     ~~Describe \code{filter_keytypes} here~~
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
function (df_pw_proteins = NULL, keytypes = c("entrezgene", "mim", 
    "uniprot", "unigene", "ensemblgene", "hprd", "hgnc"), filter_keytypes = TRUE) 
{
    require(dplyr)
    if (is.null(df_pw_proteins)) {
        stop("add_symbols_entrezids: df_pw_proteins is missing!")
    }
    df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type \%chin\% 
        "LL"] <- "entrezgene"
    df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type \%chin\% 
        "EntrezGene"] <- "entrezgene"
    df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type \%chin\% 
        "Entrez Gene"] <- "entrezgene"
    df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type \%chin\% 
        "UniProt"] <- "uniprot"
    df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type \%chin\% 
        "UniGene"] <- "unigene"
    df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type \%chin\% 
        "Ensembl"] <- "ensemblgene"
    df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type \%chin\% 
        "OMIM"] <- "mim"
    df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type \%chin\% 
        "HPRD"] <- "hprd"
    df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type \%chin\% 
        "HGNC"] <- "hgnc"
    df_pw_proteins$biopax.Gene.Symbol <- as.character(df_pw_proteins$biopax.Gene.Symbol)
    df_pw_proteins$ENTREZID <- as.integer(df_pw_proteins$ENTREZID)
    df_pw_proteins <- df_pw_proteins \%>\% expand.df.via.split.col(colToSplit = "biopax.Gene.ID")
    for (KEYTYPE in keytypes) {
        df_pw_proteins <- df_pw_proteins \%>\% add.symbols.entrezids.mygene(KEYTYPE)
    }
    df_pw_proteins <- df_pw_proteins \%>\% expand.df.via.split.col(colToSplit = c("ENTREZID", 
        "biopax.Gene.Symbol"))
    if (filter_keytypes) {
        df_pw_proteins <- df_pw_proteins \%>\% filter(biopax.Gene.ID.Type \%chin\% 
            keytypes)
    }
    return(df_pw_proteins)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

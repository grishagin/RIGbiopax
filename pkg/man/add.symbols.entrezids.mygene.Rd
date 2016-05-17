\name{add.symbols.entrezids.mygene}
\alias{add.symbols.entrezids.mygene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
add.symbols.entrezids.mygene(df_pw_proteins = NULL, KEYTYPE = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df_pw_proteins}{
%%     ~~Describe \code{df_pw_proteins} here~~
}
  \item{KEYTYPE}{
%%     ~~Describe \code{KEYTYPE} here~~
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
function (df_pw_proteins = NULL, KEYTYPE = NULL) 
{
    require(dplyr)
    require(mygene)
    if (is.null(df_pw_proteins) | is.null(KEYTYPE)) {
        stop("add_symbols_entrezids: df_pw_proteins or KEYTYPE is missing!")
    }
    KEYTYPE_rows <- which(df_pw_proteins$biopax.Gene.ID.Type \%chin\% 
        KEYTYPE)
    if (length(KEYTYPE_rows) > 0) {
        queryInput <- df_pw_proteins$biopax.Gene.ID[KEYTYPE_rows]
        queryResults <- queryMany(unique(queryInput), scopes = KEYTYPE, 
            fields = c("entrezgene", "symbol"), return.as = "DataFrame", 
            species = "human") \%>\% data.frame
        if (!("entrezgene" \%chin\% colnames(queryResults))) {
            message("Unsuccessful query for ", KEYTYPE, ".")
            return(df_pw_proteins)
        }
        queryResults <- queryResults \%>\% shrink.df.via.merge.col(colKey = "query", 
            colToMerge = c("entrezgene", "symbol"), patternToMerge = ",")
        df_pw_proteins[KEYTYPE_rows, c("ENTREZID", "biopax.Gene.Symbol")] <- queryResults[match(queryInput, 
            queryResults$query), c("entrezgene", "symbol")]
    }
    if (length(grep("netpath", unique(df_pw_proteins$biopax.Pathway.ID), 
        ignore.case = TRUE)) > 0) {
        KEYTYPE_NA_rows <- which(df_pw_proteins$biopax.Gene.ID.Type \%chin\% 
            KEYTYPE & is.na(df_pw_proteins$biopax.Gene.ID))
        if (length(KEYTYPE_NA_rows) > 0) {
            symbols <- df_pw_proteins$biopax.Component.ID[KEYTYPE_NA_rows] \%>\% 
                strsplit(split = "_") \%>\% sapply("[[", 2)
            queryResults <- queryMany(unique(symbols), scopes = "symbol", 
                fields = "entrezgene", return.as = "DataFrame", 
                species = "human") \%>\% data.frame
            if (!("entrezgene" \%chin\% colnames(queryResults))) {
                message("Unsuccessful query for symbols: ", unique(symbols))
                return(df_pw_proteins)
            }
            queryResults <- queryResults \%>\% shrink.df.via.merge.col(colKey = "query", 
                colToMerge = "entrezgene", patternToMerge = ",")
            df_pw_proteins[KEYTYPE_NA_rows, c("ENTREZID", "biopax.Gene.Symbol")] <- queryResults[match(symbols, 
                queryResults$query), c("entrezgene", "query")]
        }
    }
    return(df_pw_proteins)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

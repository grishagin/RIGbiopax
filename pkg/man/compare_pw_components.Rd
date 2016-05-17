\name{compare_pw_components}
\alias{compare_pw_components}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compare_pw_components(biopax_prot = NULL, toxdb_prot = NULL, dfrow = NULL, output = c("all", "match", "biopax", "toxdb"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{biopax_prot}{
%%     ~~Describe \code{biopax_prot} here~~
}
  \item{toxdb_prot}{
%%     ~~Describe \code{toxdb_prot} here~~
}
  \item{dfrow}{
%%     ~~Describe \code{dfrow} here~~
}
  \item{output}{
%%     ~~Describe \code{output} here~~
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
function (biopax_prot = NULL, toxdb_prot = NULL, dfrow = NULL, 
    output = c("all", "match", "biopax", "toxdb")) 
{
    dummy_df <- adjust_columns(df_to_adjust = dfrow, df_template = biopax_prot)
    pw_toxdb_id = dfrow$toxdb.Pathway.ID
    pw_biopax_id = dfrow$biopax.Pathway.ID
    if (is.na(pw_biopax_id) & is.na(pw_toxdb_id)) {
        message("compare_pw_components: both pw_biopax_id and pw_toxdb_id are NA!")
        return(dummy_df)
    }
    if (is.na(pw_biopax_id)) {
        toxdb_prot_sub <- toxdb_prot \%>\% filter(toxdb.Pathway.ID == 
            pw_toxdb_id)
        if (nrow(biopax_prot_sub) < 1) {
            message("compare_pw_components: biopax_prot_sub is NA toxdb_prot_sub has 0 rows!")
            return(dummy_df)
        }
        return(toxdb_prot_sub)
    }
    else if (is.na(pw_toxdb_id)) {
        biopax_prot_sub <- biopax_prot \%>\% filter(biopax.Pathway.ID == 
            pw_biopax_id)
        if (nrow(biopax_prot_sub) < 1) {
            message("compare_pw_components: pw_toxdb_id is NA and biopax_prot_sub has 0 rows!")
            return(dummy_df)
        }
        return(biopax_prot_sub)
    }
    biopax_prot_sub <- biopax_prot \%>\% filter(biopax.Pathway.ID == 
        pw_biopax_id & toxdb.Pathway.ID == pw_toxdb_id)
    toxdb_prot_sub <- toxdb_prot \%>\% filter(toxdb.Pathway.ID == 
        pw_toxdb_id) \%>\% mutate(biopax.Pathway.Name = as.character(dfrow$biopax.Pathway.Name), 
        biopax.Pathway.ID = as.character(dfrow$biopax.Pathway.ID), 
        pathway.Match.Status = as.character(dfrow$pathway.Match.Status))
    if (nrow(biopax_prot_sub) < 1) {
        message("compare_pw_components: biopax_prot_sub has 0 rows!")
        return(toxdb_prot_sub)
    }
    else if (nrow(toxdb_prot_sub) < 1) {
        message("compare_pw_components: toxdb_prot_sub has 0 rows!")
        return(biopax_prot_sub)
    }
    else if (nrow(biopax_prot_sub) < 1 & nrow(toxdb_prot_sub) < 
        1) {
        message("compare_pw_components: biopax_prot_sub and toxdb_prot_sub have 0 rows!")
        return(dummy_df)
    }
    pres_toxdb.pres_biopax <- biopax_prot_sub \%>\% .[(.$ENTREZID \%in\% 
        toxdb_prot_sub$toxdb.Gene.ID), ]
    pres_toxdb.pres_biopax$toxdb.Gene.ID <- toxdb_prot_sub$toxdb.Gene.ID[match(pres_toxdb.pres_biopax$ENTREZID, 
        toxdb_prot_sub$toxdb.Gene.ID)]
    pres_toxdb.pres_biopax$toxdb.Gene.Symbol <- toxdb_prot_sub$toxdb.Gene.Symbol[match(pres_toxdb.pres_biopax$ENTREZID, 
        toxdb_prot_sub$toxdb.Gene.ID)]
    pres_biopax.miss_toxdb <- biopax_prot_sub \%>\% .[!(.$ENTREZID \%in\% 
        toxdb_prot_sub$toxdb.Gene.ID), ]
    pres_toxdb.miss_biopax <- toxdb_prot_sub \%>\% .[!(.$toxdb.Gene.ID \%in\% 
        biopax_prot_sub$ENTREZID), ]
    if (output == "all") {
        to_return <- rbind(pres_toxdb.pres_biopax, pres_biopax.miss_toxdb, 
            pres_toxdb.miss_biopax)
        return(to_return)
    }
    else if (output == "match") {
        return(pres_toxdb.pres_biopax)
    }
    else if (output == "biopax") {
        return(pres_biopax.miss_toxdb)
    }
    else if (output == "toxdb") {
        return(pres_toxdb.miss_biopax)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{add.MULT.symbols.entrezids}
\alias{add.MULT.symbols.entrezids}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Add Symbols and Entrez IDs to List of Genes
}
\description{
Internal function. Add symbols and entrez IDs to a list of genes.
}
\usage{
add.MULT.symbols.entrezids(df_pw_proteins = NULL
                           ,keytypes = c("entrezgene"
                                         , "mim"
                                         , "uniprot"
                                         , "unigene"
                                         , "ensemblgene"
                                         , "hprd"
                                         , "hgnc")
                           ,filter_keytypes = TRUE)
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
Ivan Grishagin
}
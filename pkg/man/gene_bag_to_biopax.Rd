\name{gene_bag_to_biopax}
\alias{gene_bag_to_biopax}
\title{
Convert a "Gene Bag" Data Table to BioPAX Object
}
\description{
Converts a data table with pathway and gene IDs (a "gene bag") to a BioPAX object. Has an option of recording the BioPAX object to file.
}
\usage{
gene_bag_to_biopax(gene_df
                   ,filename=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_df}{
  data table with mandatory columns: \code{Pathway.Name}, \code{Pathway.ID}, \code{Gene.Symbol}, \code{Gene.ID}.
  }
  \item{filename}{
  output filename. If \code{NULL} (default), does not output to file.
  }
   
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ivan Grishagin
}
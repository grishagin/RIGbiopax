\name{merge_gene_comparison_results}
\alias{merge_gene_comparison_results}
\title{
Merge Multiple Gene Comparison Results into One File
}
\description{
Merges multiple Excel files containing gene comparison results into one Excel file. 
}
\usage{
merge_gene_comparison_results(dir=NULL
                              ,pattern=NULL
                              ,filename=NULL)
}
\arguments{
  \item{dir}{
  a directory with \code{*.xlsx} files.
}
  \item{pattern}{
  filename pattern to search for. All files whose names satisfy have the pattern will be merged.
}
  \item{filename}{
  name for the output file.
}
}

\author{
Ivan Grishagin
}


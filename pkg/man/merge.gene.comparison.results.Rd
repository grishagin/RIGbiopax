\name{merge.gene.comparison.results}
\alias{merge.gene.comparison.results}
\title{
Merge Multiple Gene Comparison Results into One File
}
\description{
Merge multiple Excel files containing gene comparison results into one BioPAX file.
}
\usage{
merge.gene.comparison.results(dir=NULL
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


\name{merge.biopax.files}
\alias{merge.biopax.files}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Merge Multiple BioPAX Files into One File
}
\description{
Merge multiple BioPAX files into one BioPAX file, while preserving all unique relationships between pathway components.
}
\usage{
merge.biopax.files(source_dir=NULL
                   ,change_ids=TRUE
                   ,write_to_file=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{source_dir}{
  a directory with \code{*.owl} files.
}
  \item{change_ids}{
  logical. Add file identifiers to every id in that file? Useful, if all files employ the same id system.
}
  \item{write_to_file}{
  logical. Write BioPAX object to file? CAUTION: takes a long time (likely, a few hours).
}
}
\details{
Returns a BioPAX object.
}

\author{
Ivan Grishagin
}


\name{MAIN_compare_toxdb_biopax}
\alias{MAIN_compare_toxdb_biopax}
\title{
Compare Proteins Extracted from BioPAX to ToxDB Proteins
}
\description{
Main function of the package. Extract all proteins ("genes") from a BioPAX source, convert their IDs to gene IDs and Symbols, and then compare those to a supplied gene list.
}
\usage{
compare_toxdb_biopax(work_dir=NULL
                     ,source_owl_dir=NULL
                     ,pw_matchup_file="./_source_files/pathways_matched_to_sources_v010.xlsx"
                     ,toxdb_genes_file="./_source_files/toxdb_pathways_15Jun_edIG_2016-08-04.txt"
                     ,source_name=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{work_dir}{
  directory, where all source files/folders are located and will be written to.
}
  \item{source_owl_dir}{
  directory with \code{*.owl} files. If no directory provided, a prompt will pop up. See \code{load.biopax} function for more information.
}
  \item{pw_matchup_file}{
  file that provides match-up between toxdb pathway names and BioPAX pathway names. 
}
\item{toxdb_genes_file}{
  file that has a table of all genes listed for all pathways. 
}
  \item{source_name}{
  name of the BioPAX source.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}

\author{
Ivan Grishagin
}

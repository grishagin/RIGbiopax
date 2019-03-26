\name{genes_and_interactions_to_biopax}
\alias{genes_and_interactions_to_biopax}
\title{
Convert a Data Table with Genes and Interactions to BioPAX Object
}
\description{
Converts a data table with pathway and gene IDs and interactions between them to a BioPAX object. 
}
\usage{
genes_and_interactions_to_biopax(gene_df
                               ,filename = NULL
                               ,cols = list(pwid = "toxdb.Pathway.ID"
                                        ,srcpwname = "source.Pathway.Name"
                                        ,first_dispname = "first_prot_sym"
                                        ,first_name = "first_name"
                                        ,first_genesym = "first_gene_sym"
                                        ,first_geneid = "first_gene_id"
                                        ,second_dispname = "second_prot_sym"	
                                        ,second_name = "second_name"
                                        ,second_genesym = "second_gene_sym"
                                        ,second_geneid = "second_gene_id"
                                        ,pmid_author = "pmid_author"
                                        ,rxn_class = "reaction_class"
                                        ,ctrl_class = "control_class"
                                        ,ctrl_type = "control_type"
                                        ,second_left_right = "second_left_right"
                                        )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_df}{
  data table with certain mandatory columns (see \code{cols} argument). Briefly, takes in a table with first and second gene parameters (gene ids, names, etc.) and reaction information (control class and type, reaction class, etc.).
  }
  \item{filename}{
  output filename. If \code{NULL} (default), does not output to file.
  }
   \item{cols}{
  \code{gene_df} columns assigned to self-explanatory variables used within the function.
  }
   
}
\details{
Has an option of recording the BioPAX object to file. All necessary columns are defined in \code{cols} variable. Returns a \code{biopax} object.
}
\value{
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ivan Grishagin
}
\name{compare_toxdb_biopax}
\alias{compare_toxdb_biopax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compare_toxdb_biopax(source_dir = "D:/Dropbox/Rancho/NCATS/ToxDB/_owls", pw_matchup_file = "./_source_files/pathways_matched_to_sources_v008.xlsx", source_name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{source_dir}{
%%     ~~Describe \code{source_dir} here~~
}
  \item{pw_matchup_file}{
%%     ~~Describe \code{pw_matchup_file} here~~
}
  \item{source_name}{
%%     ~~Describe \code{source_name} here~~
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
function (source_dir = "D:/Dropbox/Rancho/NCATS/ToxDB/_owls", 
    pw_matchup_file = "./_source_files/pathways_matched_to_sources_v008.xlsx", 
    source_name = NULL) 
{
    biopax_source_names <- c("BioCarta", "KEGG", "NCI-Nature", 
        "NetPath", "Reactome", "Science Signaling", "Wiki Pathways")
    require(RIGessentials)
    pkg <- c("plyr", "dplyr", "rBiopaxParser", "mygene", "readxl", 
        "openxlsx")
    loadPackages(pkg)
    prepareSession("D:/Dropbox/Rancho/NCATS/ToxDB/", nolocale = FALSE)
    if (!(source_name \%in\% biopax_source_names)) {
        source_name <- tkradio_from_vect(biopax_source_names, 
            "Select BioPAX Source.")
    }
    pathways_per_source <- read_excel(path = pw_matchup_file, 
        col_types = rep("text", 10), sheet = 1) \%>\% dplyr::filter(!is.na(biopax.Pathway.Name)) \%>\% 
        dplyr::select(toxdb.Pathway.ID, toxdb.Pathway.Name, biopax.Pathway.Name, 
            pathway.Match.Status, Source) \%>\% arrange(Source, 
        toxdb.Pathway.ID)
    all_pathways <- pathways_per_source \%>\% filter(Source == 
        source_name)
    all_pathways$toxdb.Pathway.Name <- tolower(all_pathways$toxdb.Pathway.Name)
    toxdb <- load.toxdb.genes.per.source(source_name = source_name, 
        source_dir = source_dir, all_pathways = all_pathways)
    owl_biopax <- load.biopax(source_name = source_name, source_dir = source_dir)
    pw_biopax <- load.biopax.pathways(owl_biopax = owl_biopax)
    pw_biopax$biopax.Pathway.Name <- pw_biopax$biopax.Pathway.Name \%>\% 
        trimws
    mismatch_pathways <- all_pathways[!(all_pathways$biopax.Pathway.Name \%in\% 
        pw_biopax$biopax.Pathway.Name) | !(all_pathways$toxdb.Pathway.Name \%in\% 
        toxdb$toxdb.Pathway.Name), ]
    if (nrow(mismatch_pathways) > 0) {
        message("Some pathways from the list do not match those in biopax or toxdb gene list!")
        print(mismatch_pathways)
        View(mismatch_pathways)
        stop("Aborting.")
    }
    biopax_pw_missing_from_list <- pw_biopax[!(pw_biopax$biopax.Pathway.Name \%chin\% 
        all_pathways$biopax.Pathway.Name), ] \%>\% mutate(pathway.Match.Status = "Missing", 
        toxdb.Pathway.ID = NA, toxdb.Pathway.Name = NA, Source = source_name)
    all_pathways_pwid <- add.biopax.ids(all_pathways = all_pathways, 
        pw_biopax = pw_biopax)
    all_pathways_pwid <- rbind.data.frame(all_pathways_pwid, 
        biopax_pw_missing_from_list)
    df_pw_proteins <- adply(.data = all_pathways_pwid, .margins = 1, 
        .fun = function(x) {
            add_db_ids(owl_biopax = owl_biopax, pw_id = x$biopax.Pathway.ID)
        }) \%>\% mutate(biopax.Gene.Symbol = NA, ENTREZID = NA) \%>\% 
        data.frame
    df_pw_proteins_annot <- add.MULT.symbols.entrezids(df_pw_proteins = df_pw_proteins, 
        filter_keytypes = TRUE)
    df_pw_proteins_annot <- adjust_columns(df_to_adjust = df_pw_proteins_annot, 
        df_template = toxdb)
    toxdb <- adjust_columns(df_to_adjust = toxdb, df_template = df_pw_proteins_annot)
    toxdb$pathway.Match.Status <- all_pathways_pwid$pathway.Match.Status[match(toxdb$toxdb.Pathway.ID, 
        all_pathways_pwid$toxdb.Pathway.ID)]
    comparison_results <- adply(.data = as.data.frame(all_pathways_pwid), 
        .margins = 1, .fun = function(x) {
            compare_pw_components(biopax_prot = df_pw_proteins_annot, 
                toxdb_prot = toxdb, dfrow = x, output = "all")
        })
    openxlsx:::write.xlsx(all_pathways_pwid, file = paste(Sys.Date(), 
        "pathways", source_name, "toxdb_biopax.xlsx", sep = "_"), 
        col.names = TRUE, row.names = FALSE)
    openxlsx:::write.xlsx(comparison_results, file = paste(Sys.Date(), 
        "genes", source_name, "toxdb_biopax.xlsx", sep = "_"), 
        col.names = TRUE, row.names = FALSE, keepNA = TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

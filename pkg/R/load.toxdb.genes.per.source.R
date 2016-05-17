load.toxdb.genes.per.source <-
function(source_name
             ,source_dir=getwd()
             ,toxdb_genes_file="toxdb_pathways_15Jun.txt"
             ,all_pathways){
        #list of genes and pathways
        #and select only pertaining columns
        options(stringsAsFactors = FALSE)
        
        toxdb<-read.table(file.path(source_dir
                                    ,toxdb_genes_file)
                          ,header = TRUE
                          ,quote=""
                          ,sep="\t") %>% 
            dplyr::select(toxdb.Pathway.ID=Pathway.ID
                          ,toxdb.Pathway.Name=Pathway.Name
                          ,toxdb.Gene.ID=Gene.ID
                          ,toxdb.Gene.Symbol=Gene.Symbol) 
        #pathway name and gene symbol identifiers -> to lower case
        toxdb$toxdb.Pathway.Name<-
            tolower(toxdb$toxdb.Pathway.Name)
        
        #get a subset of toxdb for just one source (as featured in Ruili's list of curated pathways)
        toxdb<-
            toxdb %>%
            subset(.$toxdb.Pathway.ID %in% all_pathways$toxdb.Pathway.ID) %>%
            mutate(Source=source_name)
        
        return(toxdb)
    }

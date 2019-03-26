load_inxight_genes_per_source<-
    function(source_name
             ,toxdb_genes_file=NULL
             ,all_pathways){
        #' @title
        #' Load Dataframe with Genes/Pathways
        #' @description 
        #' Load a table with all genes listed for each of the pathways.
        #' @param source_name Name of the BioPAX source.
        #' @param toxdb_genes_file File that has a table of all genes listed for all pathways. 
        #' @param all_pathways Dataframe with all desired pathways.
        
        #' @author 
        #' Ivan Grishagin
        
        #list of genes and pathways
        #and select only pertaining columns
        options(stringsAsFactors = FALSE)
        
        toxdb<-read.table(toxdb_genes_file
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

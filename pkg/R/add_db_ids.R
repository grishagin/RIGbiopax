add_db_ids <-
function(owl_biopax=NULL
             ,pw_id=NULL
             ,biopaxlevel=3){
        
        #function obtains a dataframe of db names and ids 
        #for required pathway
        if(is.null(owl_biopax) || 
           is.null(pw_id)
        ){
            stop("Some of the required parameters were not provided.")
        } else if (is.na(pw_id)){
            results_df<-
                data.frame(biopax.Component.ID=NA
                           ,biopax.Gene.ID.Type=NA
                           ,biopax.Gene.ID=NA)
            return(results_df)
        }
        
        
        require(tidyr)
        require(dplyr)
        require(data.table)
        require(rBiopaxParser)
        
        #component ids
        pw_components_ids<-
            listPathwayComponents(biopax=owl_biopax
                                  ,id=pw_id
                                  ,returnIDonly=TRUE
                                  ,biopaxlevel=biopaxlevel)
        
        #get component instances
        pw_components<-
            selectInstances(biopax=owl_biopax
                            ,id=pw_components_ids
                            ,includeReferencedInstances=TRUE)
        
        #determine Xref class name
        print(pw_id)
        
        #get a dataframe of db names and db ids
        #also add pw name and pw id columns
        results_df<-
            pw_components %>%
            filter(class=="RelationshipXref" | class=="UnificationXref") %>%
            dplyr::select(biopax.Component.ID=id
                          ,property
                          ,property_value) %>%
            spread(property
                   ,property_value)
        
        if(nrow(results_df)<1){
            return(data.frame(biopax.Component.ID=NA
                              ,biopax.Gene.ID.Type=NA
                              ,biopax.Gene.ID=NA))
        }
        
        setnames(results_df,"db","biopax.Gene.ID.Type")
        setnames(results_df,"id","biopax.Gene.ID")
        
        #replace empty values with NA
        results_df$biopax.Gene.ID[results_df$biopax.Gene.ID %chin% ""]<-
            NA
        
        results_df<-
            results_df %>%
            dplyr::select(biopax.Component.ID
                          ,biopax.Gene.ID.Type
                          ,biopax.Gene.ID
            )
        
        return(results_df)
    }

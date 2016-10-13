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
        
        #now, some biopax sources are prepared from multiple files
        #which sometimes do not have explicit relationships between components
        #and genes/proteins (i.e. these pathways are simply "gene bags")
        #to include them, first determine if there's a pathway reference
        
        #deprecated
        # if(length(grep("_pwref_",pw_id))>0){
        #     pwref<-
        #         strsplit(pw_id,"_pwref_") %>%
        #         unlist %>%
        #         .[length(.)]
        #     #get all instances with said pwref
        #     pw_components_ids<-
        #         owl_biopax$dt$id %>%
        #         grep(pwref
        #              ,.
        #              ,value = TRUE)
        #     
        # } else {
            #if no such reference has been found
            #look for component ids
        pw_components_ids<-
            listPathwayComponents(biopax=owl_biopax
                                  ,id=pw_id
                                  ,returnIDonly=TRUE
                                  ,biopaxlevel=biopaxlevel)
        #}
        if (is.null(pw_components_ids) | 
            is.na(pw_components_ids)){
            stop("add_db_ids: could not find any components in pathway "
                 ,pw_id
                 ,"!")
        }
    
        #get component instances
        pw_components<-
            selectInstances(biopax=owl_biopax
                            ,id=pw_components_ids
                            ,includeReferencedInstances=TRUE) %>%
            unique
            
        
        #determine Xref class name
        #print(pw_id)
        
        #get a dataframe of db names and db ids
        #also add pw name and pw id columns
        results_df<-
            pw_components %>%
            filter(class %in% "RelationshipXref" | class %in% "UnificationXref") %>%
            dplyr::select(biopax.Component.ID=id
                          ,property
                          ,property_value)
        
        if(nrow(results_df)<1){
            return(data.frame(biopax.Component.ID=NA
                              ,biopax.Gene.ID.Type=NA
                              ,biopax.Gene.ID=NA))
        } else {
            results_df<-
                results_df %>%
                spread(property
                       ,property_value)
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

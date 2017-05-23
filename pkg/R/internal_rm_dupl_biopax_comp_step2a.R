internal_rm_dupl_biopax_comp_step2a<-
    function(biopax){
        #this function sorts out all duplicate instances
        #referring to the term-db-id, etc.
        #such as xref, term, evidenceCode, featureLocation, etc.
        #needs to be applied iteratively to cover all tiers of nesting
        
        #if a components with different ids have only certain properties (keep_props) and 
        #the other fields of these records are the same (i.e. pavs and/or property values) -- they are the same
        orig_biopax<-
            copy(biopax)
        keep_props<-
            c("term"
              ,"xref"
              ,"evidenceCode"
              ,"featureLocation"
              ,"modificationType"
              ,"sequenceIntervalBegin"
              ,"sequenceIntervalEnd")
        
        #find all ids with duplicated pavs (which are refs to other components)
        dupl_pav_ids<-
            biopax$dt[property_attr_value %in% property_attr_value[duplicated(property_attr_value)] &
                          property_attr_value %in% id]$id
        
        #of those, take those, which have ONLY allowed properties
        ids_to_proc<-
            biopax$dt[id %in% dupl_pav_ids
                      ,.(to_take=all(property %in% keep_props))
                      ,by=id][,.(id=id[to_take])]$id
        
        #if there are duplicate ids, process
        if(length(ids_to_proc)>1){
            #ids with their replacement values
            ids_to_proc_df<-
                biopax$dt[id %in% ids_to_proc &
                              property_attr_value %in% id
                          ,.(from=id
                             ,to=id[1])
                          ,by=property_attr_value][from!=to
                                                   ,.(from
                                                      ,to)] %>% 
                unique %>% 
                #some properties (e.g. Begin and End) have the same ids
                #so some from ids will be the same -- exclude those
                .[!duplicated(from)]
            
            
            #now replace all those ids in the dataframe
            #in id and pav columns
            biopax$dt$id<-
                biopax$dt$id %>% 
                mapvalues(from=ids_to_proc_df$from
                          ,to=ids_to_proc_df$to)
            biopax$dt$property_attr_value<-
                biopax$dt$property_attr_value %>% 
                mapvalues(from=ids_to_proc_df$from
                          ,to=ids_to_proc_df$to
                          #,warn_missing = FALSE
                          )
            
            #take only unique rows
            biopax$dt<-
                biopax$dt %>% 
                unique
        } else {
            biopax<-
                orig_biopax
        }
        #clean-up
        rm(orig_biopax)
        
        return(biopax)
        
    }
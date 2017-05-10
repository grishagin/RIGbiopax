internal_rm_dupl_biopax_comp_step2<-
    function(biopax){
        #if a components with different ids have only certain properties (keep_props) and 
        #the other fields of these records are the same (i.e. pavs and/or property values) -- they are the same
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
        
        #if no ids, return original biopax
        if(length(ids_to_proc)<1){
            #NB -- at this point original biopax != this biopax
            #as those operations actually changed indexing in the data table
            #but as far as data is concerned, they are the same
            return(biopax)
        }

        #ids with their replacement values
        ids_to_proc_df<-
            biopax$dt[id %in% ids_to_proc &
                          property_attr_value %in% id
                      ,.(from=id
                         ,to=id[1])
                      ,by=property_attr_value][from!=to] %>% 
            unique
        
        #now replace all those ids in the dataframe
        #in id and pav columns
        biopax$dt$id<-
            biopax$dt$id %>% 
            mapvalues(from=ids_to_proc_df$from
                      ,to=ids_to_proc_df$to)
        biopax$dt$property_attr_value<-
            biopax$dt$property_attr_value %>% 
            mapvalues(from=ids_to_proc_df$from
                      ,to=ids_to_proc_df$to)
        
        #take only unique rows
        biopax$dt<-
            biopax$dt %>% 
            unique
        
        #try the same function again recursively
        #until it returns the same result
        biopax2<-NULL
        #will run at least once
        while(!identical(biopax,biopax2)){
            #store old biopax in a new variable
            biopax2<-
                copy(biopax)
            #try function on current biopax
            biopax<-
                biopax %>% 
                internal_rm_dupl_biopax_comp_step2
            #if any other mods are avaialble, 
            #it will be different from the old one
        }
        #clean up
        rm(biopax2)
        
        return(biopax)
        
    }
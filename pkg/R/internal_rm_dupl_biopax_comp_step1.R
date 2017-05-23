internal_rm_dupl_biopax_comp_step1<-
    function(biopax){
        #this function sorts out all duplicate term-db-id, db-id, and position-status values
        
        #create an explicit copy, such that the original object is not affected
        biopax_dt<-
            copy(biopax$dt)
        
        keep_props<-
            c("term"
              ,"db"
              ,"id"
              ,"positionStatus"
              ,"sequencePosition")
        
        #merge class, property, and property value
        #only for instances with "good" properties
        #i.e. explicitly exclude all other properties with same id
        biopax_dt[!id %in% id[!property %in% keep_props]
                  ,combo1 := paste0(class
                                    ,property
                                    ,property_value)]
        biopax_dt<-
            biopax_dt[order(id
                            ,combo1)]
        
        #merge all such values within one id
        biopax_dt[combo1!="NA"
                  ,combo2 := paste0(combo1
                                    ,collapse="")
                  ,by=id]
        
        #now for each supercombo (combo2)
        #assign replacement ids -- just pick the first of the group of ids
        #with the same combo2
        allid_by_combo2<-
            biopax_dt[!is.na(combo2)
                      ,.(from=id
                         ,to=id[1])
                      ,by=combo2][from!=to] %>% 
            unique
        
        #now replace all those ids in the dataframe
        #in id and pav columns
        biopax_dt$id<-
            biopax_dt$id %>% 
            mapvalues(from=allid_by_combo2$from
                      ,to=allid_by_combo2$to)
        biopax_dt$property_attr_value<-
            biopax_dt$property_attr_value %>% 
            mapvalues(from=allid_by_combo2$from
                      ,to=allid_by_combo2$to
                      ,warn_missing = FALSE)
        
        #take only unique rows
        biopax_dt<-
            biopax_dt %>% 
            unique
        
        #remove auxiliary columns
        biopax_dt$combo1<-NULL
        biopax_dt$combo2<-NULL
        
        #replace original datatable with the modified one
        biopax$dt<-
            biopax_dt
        #clean up
        rm(biopax_dt)
        
        return(biopax)
        
    }
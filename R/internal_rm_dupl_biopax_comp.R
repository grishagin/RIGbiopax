internal_rm_dupl_biopax_comp<-
  function(biopax
           ,logical_subset_vect){
    
    #' @keywords internal
    
    #make a copy of the biopax dt
    biopax_dt<-
      copy(biopax$dt)
    
    #merge class, property, property_attr_value, and property_value
    #only for desired instances 
    biopax_dt[logical_subset_vect
              ,combo1 := paste0(class
                                ,property
                                ,property_attr_value
                                ,property_value)]
    
    #order by id, then combo1 -- for consistency later
    biopax_dt<-
      biopax_dt[order(id
                      ,combo1)]
    
    #merge all non-NA combo1 values within one id
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
    # print(allid_by_combo2[from %in% to])
    # print(allid_by_combo2[to %in% from])
    # print(allid_by_combo2)
    
    #if no such rows found -- return original biopax
    if(nrow(allid_by_combo2)<1){
      rm(biopax_dt)
      return(biopax)
    }
    
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
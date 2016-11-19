unify_biopax_ids<-
    function(biopax_dt
             ,idtag=NULL
             ,exclude_id_pattern="bioplanet"
             ,exclude_class="Pathway"){
        #standardizes representation of all biopax ids
        #take unique classes with ids combinations
        #doesn't touch ids with exclude_id_pattern of class exclude_class
        st<-Sys.time()
        class_id<-
            biopax_dt %>%
            dplyr::select(class,id) %>%
            unique 
        if(!is.null(exclude_id_pattern)){
            class_id<-
                class_id  %>%
                #exclude instances with ids featuring a given pattern
                #AND belonging to the class to be excluded
                filter(!(grepl(exclude_id_pattern
                               ,id) &
                             class %in% exclude_class)) 
        }
       
        
        if (nrow(class_id)<1){
            message("unify_biopax_ids: no ids to replace!")
            return(biopax_dt)
        }
        #some entities refer to instances of multiple different classes
        #in that case, we'll replace class with "PhysicalEntity"
        dupl_ids<-
            class_id$id[duplicated(class_id$id)]
        class_id$class[class_id$id %in% dupl_ids]<-
            "PhysicalEntity"
        
        #make new ids by merging class with its "order number"
        class_id<-
            class_id %>%
            unique %>%
            mutate(newid = paste(idtag
                                 ,class
                                 ,RIGbiopax:::internal_seq_along_find_reps(class)
                                 ,sep=""))
        #replace ids
        biopax_dt$id<-
            biopax_dt$id %>%
            mapvalues(from = class_id$id
                      ,to = class_id$newid)
        #remove ids that are not referenced in the original dataframe
        class_id<-
            class_id %>%
            filter(id %in% biopax_dt$property_attr_value)

        #replace prop attr values (references)
        biopax_dt$property_attr_value<-
            biopax_dt$property_attr_value %>%
            mapvalues(from = class_id$id
                      ,to = class_id$newid)
        et<-Sys.time()
        #print(et-st)
        
        
        return(biopax_dt)
        
    }
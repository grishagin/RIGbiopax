unify_biopax_ids<-
    function(biopax_dt
             ,idtag=NULL){
        #standardizes representation of all biopax ids
        #take unique classes with ids combinations
        st<-Sys.time()
        class_id<-
            biopax_dt %>%
            dplyr::select(class,id) %>%
            filter(class!="Pathway") %>%
            unique
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
        #replace prop attr values (references)
        biopax_dt$property_attr_value<-
            biopax_dt$property_attr_value %>%
            mapvalues(from = class_id$id
                      ,to = class_id$newid)
        et<-Sys.time()
        print(et-st)
        
        
        return(biopax_dt)
        
    }
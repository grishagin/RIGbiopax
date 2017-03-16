unify_biopax_ids<-
    function(biopax_dt
             ,idtag=NULL
             ,exclude_id_pattern="inxight_pathways"
             ,exclude_class="Pathway"){
        
        #' @title
        #' Uniformly Format IDs in a BioPAX-Style Data Table
        #' @description 
        #' Converts IDs in a BioPAX-style data table to a uniform format.
        #' @details 
        #' ID unification is conducted by appending a class of a component to its throughout number.
        #' E.g. for a component of class Protein, #11 from the top of the table, ID will be "Protein11".
        #' @param biopax_dt BioPAX-style data table.
        #' @param idtag String to append to each id (defaults to \code{NULL}).
        #' @param exclude_id_pattern Exclude components with such pattern in their IDs and of class \code{exclude_class} from ID unification.
        #' @param exclude_class Exclude components with such class and \code{exclude_id_pattern} from ID unification.
        
        #' @author 
        #' Ivan Grishagin
        
        #standardizes representation of all biopax ids
        #take unique classes with ids combinations
        #doesn't touch ids with exclude_id_pattern of class exclude_class
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
                               ,id) 
                         & class %in% exclude_class)) 
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
                                 ,internal_seq_along_find_reps(class)
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

        return(biopax_dt)
        
    }
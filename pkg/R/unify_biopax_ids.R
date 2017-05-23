unify_biopax_ids<-
    function(biopax
             ,idtag=NULL
             ,exclude_id_pattern="inxight_pathways"
             ,exclude_class="Pathway"){
        
        #' @title
        #' Uniformly Format IDs in BioPAX
        #' @description 
        #' Converts IDs in a BioPAX object or a BioPAX-style data table to a uniform format.
        #' @details 
        #' ID unification is conducted by appending a class of a component to its throughout number.
        #' E.g. for a component of class Protein, #11 from the top of the table, ID will be "Protein11". 
        #' In addition, component classes are corrected as well to ensure that components with the same id have the same class.
        #' @param biopax BioPAX object or a BioPAX-style data table.
        #' @param idtag String to append to each id (defaults to \code{NULL}).
        #' @param exclude_id_pattern Exclude components with such pattern in their IDs and of class \code{exclude_class} from ID unification.
        #' @param exclude_class Exclude components with such class and \code{exclude_id_pattern} from ID unification.
        
        #' @author 
        #' Ivan Grishagin
        
        return_biopax<-
            FALSE
        
        #check input
        if("biopax" %in% class(biopax)){
            return_biopax<-
                TRUE
            biopax_dt<-
                biopax$dt
        } else if (is.data.table(biopax)){
            biopax_dt<-
                biopax
        } else if (is.data.frame(biopax)){
            biopax_dt<-
                biopax %>%
                as.data.table
        } else {
            stop("unify_biopax_ids: something's wrong with the supplied 'biopax' argument!")
        }
        
        
        #take unique classes with ids combinations
        class_id<-
            biopax_dt[,.(class
                         ,id)] %>%
            unique 
        
        if(!is.null(exclude_id_pattern)){
            #exclude instances with ids featuring a given pattern
            #AND belonging to the class to be excluded
            class_id<-
                class_id[!(grepl(exclude_id_pattern
                                 ,id) 
                           & class %in% exclude_class)]
                
        }
       
        
        if (nrow(class_id)<1){
            message("unify_biopax_ids: no ids to replace!")
            return(biopax)
        }
        #some entities refer to instances of multiple different classes
        
        #for those dupl ids that have class ModificationFeature and FragmentFeature
        #replace them with just ModificationFeature 
        class_id[id %in% id[duplicated(id)]][id %in% id[class=="ModificationFeature"]]$class<-
            "ModificationFeature"
        class_id<-
            class_id %>%
            unique 
        
        #in leftover cases, we'll replace class with "PhysicalEntity"
        class_id[id %in% id[duplicated(id)]]$class<-
            "PhysicalEntity"
        class_id<-
            class_id %>%
            unique 
        
        class_id[,newid := paste0(idtag
                                  ,class
                                  ,1:length(id))
                 ,by=class]
        
        #replace ids
        biopax_dt$id<-
            biopax_dt$id %>%
            mapvalues(from = class_id$id
                      ,to = class_id$newid)
        
        #replace classes
        biopax_dt[id %in% class_id$newid]$class<-
            biopax_dt[id %in% class_id$newid]$id %>%
            match(class_id$newid) %>% 
            class_id$class[.]
           
        
        #remove ids that are not referenced in the original dataframe
        #from class_id dataframe
        class_id<-
            class_id[id %in% biopax_dt$property_attr_value]

        #replace prop attr values (references)
        biopax_dt$property_attr_value<-
            biopax_dt$property_attr_value %>%
            mapvalues(from = class_id$id
                      ,to = class_id$newid)
        
        #clean up
        class_id<-NULL
        
        #if biopax object was originally supplied,
        #convert the modified data table back into the biopax object
        if(return_biopax){
            biopax$dt<-
                biopax_dt
        } else {
            biopax<-
                biopax_dt
        }
        return(biopax)
        
    }
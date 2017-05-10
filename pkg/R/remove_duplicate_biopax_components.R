remove_duplicate_biopax_components<-
    function(biopax){
        #' @title 
        #' Remove Duplicate Components from BioPAX Object
        #' @description 
        #' BioPAX objects frequently have components that are identical in everything, but id. 
        #' This function removes all such duplicates, and changes all property_attr_values accordingly.
        #' @author 
        #' Ivan Grishagin
        
        biopax<-
            biopax %>% 
            #fix db/id positionStatus/sequencePosition and term entries
            internal_rm_dupl_biopax_comp_step1 %>% 
            #fix all properties that link to them
            internal_rm_dupl_biopax_comp_step2 %>% 
            #iteratively fix physical entities
            internal_rm_dupl_biopax_comp_step3
        
        return(biopax)
    }
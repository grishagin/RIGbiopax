remove_deadend_refs<-
    function(biopax){
        #' @title 
        #' Remove Dead-End References from BioPAX
        #' @description 
        #' Takes in a BioPAX file, and removes those instances that refer 
        #' to instances that do not exist.
        #' @param biopax A BioPAX object.
        #' @author 
        #' Ivan Grishagin
        
        to_remove_logi<-
            biopax$dt[,.(to_remove=
                             property_attr=="rdf:resource" & 
                             !property_attr_value %in% id)]$to_remove
        
        if(sum(to_remove_logi)>0){
            #if some rows need to be removed
            #remove them
            biopax$dt<-
                biopax$dt[!to_remove_logi]
            #and repeat to avoid situation where 
            #removal of dead references creates new ones
            biopax<-
                biopax %>% 
                remove_deadend_refs
        }

        return(biopax)
        
        }
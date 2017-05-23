internal_rm_dupl_biopax_comp_step2<-
    function(biopax){
        #this function is a wrapper for another function
        #which sorts out all duplicate instances
        #referring to the term-db-id, etc.
        #such as xref, term, evidenceCode, featureLocation, etc.
        
        #try the same function recursively
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
                internal_rm_dupl_biopax_comp_step2a
            #if any other mods are avaialble, 
            #it will be different from the old one
        }
        #clean up
        rm(biopax2)
        
        return(biopax)
        
    }
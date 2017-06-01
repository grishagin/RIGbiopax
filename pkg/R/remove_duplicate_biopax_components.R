remove_duplicate_biopax_components<-
    function(biopax){
        #' @title 
        #' Remove Duplicate Components from BioPAX Object
        #' @description 
        #' BioPAX objects frequently have components that are identical in everything, but id. 
        #' This function removes all such duplicates, and changes all property_attr_values accordingly.
        #' @details 
        #' The biopax is processed iteratively.\cr
        #' First, two sets of properties are used to identify and remove:\cr
        #' 1) duplicate term-db-id, db-id, and position-status instances;\cr
        #' 2) duplicate instances that refer to instance listed in 1) -- 
        #' they have properties xref, term, evidenceCode, featureLocation, 
        #' modificationType, sequenceIntervalBegin, sequenceIntervalEnd. \cr\cr
        #' Second, two sets of classes are used to identify and remove:\cr
        #' 1) duplicate Protein, SmallMolecule, PhysicalEntity, RNA, DNA instances;\cr
        #' 2) duplicate Complex instances.\cr
        #' 
        #' @author 
        #' Ivan Grishagin
        
        #strip hashes from the biopax property attr value column
        biopax$dt$property_attr_value<-
            biopax$dt$property_attr_value %>% 
            striphash
        
        #first, clean-up entities by properties
        #then by class
        for(keep_list_name in c("property","class")){
            
            message("Removal of duplicate biopax components: using "
                    ,toupper(keep_list_name)
                    ," started...")
            
            #list of properties or classes to preserve
            keep_list<-
                internal_keep_property_class_list(keep_list_name=keep_list_name)
            
            #for each level of the element (property or class), 
            #perform the biopax component clean-up procedure
            for (keep_list_element in keep_list){
                #auxiliary biopax to compare to
                aux_biopax<-NULL
                iteration=0
                
                while(!identical(biopax,aux_biopax)){
                    iteration<-
                        iteration+1
                    message("Removal of duplicate biopax components: iteration "
                            ,iteration
                            ," started...")
                    
                    #store old biopax in a new variable
                    aux_biopax<-
                        copy(biopax)
                    
                    #make a logical vector to subset the biopax
                    #i.e. determine all instances to be processed
                    logical_subset_vect<-
                        switch(keep_list_name
                               ,property = biopax$dt[,.(logi=!id %in% id[!property %in% keep_list_element])]$logi
                               ,class = biopax$dt[,.(logi=class %in% keep_list_element)]$logi
                        )
                    
                    #try function on "old" biopax
                    #in the next iteration it will be compared to old self
                    #and if different, the procedure will be repeated
                    biopax<-
                        biopax %>% 
                        internal_rm_dupl_biopax_comp(logical_subset_vect)
                  
                    #if any other mods are avaialble, 
                    #it will be different from the old one
                    
                    message("Removal of duplicate biopax components: iteration "
                            ,iteration
                            ," is over!")
                }
                #clean up
                rm(aux_biopax)
            }
            
            message("Removal of duplicate biopax components: using "
                    ,toupper(keep_list_name)
                    ," is over!")
        }
        
        return(biopax)
        
    }
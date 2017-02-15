extract_pathways_from_biopax<-
    function(biopax
             ,pw_ids
             ,filename=NULL){
        #'@title
        #'Extract Pathways from BioPAX
        #'@description
        #'Extract pathways from a BioPAX object, then convert them into a separate biopax object.
        #'@param biopax BioPAX object.
        #'@param pw_ids Pathway ids to extract.
        #'@param filename Filename of the new biopax.
        #'
        #'@author
        #'Ivan Grishagin
        
        #check if all ids have been found in biopax
        not_in_biopax<-
            pw_ids[!pw_ids %in% listPathways(biopax)$id]
        if(length(not_in_biopax)>0){
            warning("Some pathway ids you provided have not been found in biopax: "
                    ,paste(not_in_biopax
                           ,collapse = ", ")
                    ,".\nOnly other those ids that have been found will be processed.")
            pw_ids<-
                pw_ids[!pw_ids %in% not_in_biopax]
        }
        
        if(length(pw_ids)<1){
            warning("None of the supplied pathway ids were found in the biopax! Aborting.")
            return()
        }
        
        #get all referenced instances and convert them to biopax
        new_biopax<-
            pw_ids %>% 
            selectInstances(biopax = biopax
                            ,id = .
                            ,includeReferencedInstances=TRUE) %>% 
            biopax_from_dt(filename = filename)
        return(new_biopax)
        
    }
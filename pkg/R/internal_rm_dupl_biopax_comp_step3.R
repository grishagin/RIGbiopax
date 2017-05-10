internal_rm_dupl_biopax_comp_step3<-
    function(biopax){
        #this function calls itself iteratively 
        #first, remove copies of all physical entities except Complex
        #then do same for complexes (iteratively)
        #then reactions
        #then control types
        #going higher and higher in biopax hierarchy
      
        #list of classes to keep
        keep_class_list<-
            internal_keep_class_list()
        
        for (lvl in 1:length(keep_class_list)){
            keep_class<-
                keep_class_list[[lvl]]
            biopax<-
                biopax %>% 
                internal_rm_dupl_biopax_comp_step3a(keep_class = keep_class)
            
            if("Complex" %in% keep_class){
                #if running for complex, repeat until resultant biopax come out identical
                #since there can be complexes inside complexes
                biopax2<-NULL
                
                #will run at least once
                while(!identical(biopax,biopax2)){
                    #store old biopax in a new variable
                    biopax2<-
                        copy(biopax)
                    #try function on current biopax
                    biopax<-
                        biopax %>% 
                        internal_rm_dupl_biopax_comp_step3a(keep_class = keep_class)
                    #if any other mods are avaialble, 
                    #it will be different from the old one
                }
                #clean up
                rm(biopax2)
            }
        }

        
        return(biopax)
        
    }
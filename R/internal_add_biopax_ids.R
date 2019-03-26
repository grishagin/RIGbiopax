internal_add_biopax_ids<-
    function(all_pathways
             ,pw_biopax){
        #add biopax pathway ids
        all_pathways_pwid<-
            adply(.data = as.data.frame(all_pathways)
                  ,.margins = 1
                  ,.fun = function(dfrow){
                      biopaxID<-
                          pw_biopax$biopax.Pathway.ID[pw_biopax$biopax.Pathway.Name %chin% 
                                                          dfrow$biopax.Pathway.Name]
                      dfrow<-
                          dfrow %>%
                          .[rep(1,length(biopaxID)),] %>%
                          mutate(biopax.Pathway.ID=biopaxID)
                      
                      if(length(biopaxID)>1) {
                          message("add_biopax_ids: Found "
                                  ,length(biopaxID)
                                  ," IDs for "
                                  ,unique(dfrow$biopax.Pathway.Name))
                      }
                      return(dfrow)
                      
                  })
        return(all_pathways_pwid)
    }

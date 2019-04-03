extract_pathways_from_biopax<-
  function(biopax
           ,pw_ids
           ,filename=NULL
           ,exclude_subpw_comp_patt=NULL){
    #' @export
    #' @title
    #' Extract Pathways from BioPAX
    #' @description
    #' Extract pathways from a BioPAX object, then convert them into a separate biopax object.
    #' @param biopax BioPAX object.
    #' @param pw_ids Pathway ids to extract.
    #' @param filename Filename of the new BioPAX
    #' @param exclude_subpw_comp_patt Pattern of sub-pathways, whose components will be excluded, 
    #' i.e. the main pathway will refer to them, but the references to their components will be removed.
    #' 
    #' @author
    #' Ivan Grishagin
    
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
    
    #exclude some components
    #NOT(exclude pattern IN id AND NOT id IN desired pw ids AND property is component)
    if(!is.null(exclude_subpw_comp_patt)){
      biopax$dt<-
        biopax$dt[!(grepl(exclude_subpw_comp_patt
                          ,id) &
                      !id %in% pw_ids &
                      property=="pathwayComponent")]
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
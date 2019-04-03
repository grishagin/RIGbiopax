split_pathway_into_chunks<-
  function(biopax
           ,pw_id
           ,Nchunks=2){
    #' @export
    #' @title
    #' Split BioPAX Pathway into N Chunks
    #' @description
    #' Split a given pathway supplied in a BioPAX object into N new BioPAX objects.
    #' @param biopax BioPAX object.
    #' @param pw_id ID of a pathway to split.
    #' @param Nchunks How many new BioPAX objects to produce from the original one?
    #'
    #' @author
    #' Ivan Grishagin
    
    #check if the pathway id has been found in biopax
    if(!pw_id %in% listPathways(biopax)$id){
      warning("split_pathway_into_chunks: none of the supplied pathway ids were found in the biopax! Aborting.")
      return()
    }
    
    
    
    #get ids of pathway components
    pw_comp_ids<-
      biopax %>% 
      listPathwayComponents(id = pw_id
                            ,includeSubPathways = FALSE
                            ,returnIDonly = TRUE)
    
    ############################################# split the id vector up into Nchunks
    #safeguard agains too manu chuncks
    Nchunks<-
      min(Nchunks,length(pw_comp_ids))
    
    #get split intervals based on desired number of intervals 
    quant_splits<-
      quantile(1:length(pw_comp_ids)
               ,probs=seq(0,1
                          ,1/Nchunks)
               ,type=1)
    #change the first element to zero (see below)
    quant_splits[1]<-0
    #define empty list to hold pw comp ids
    pw_comp_list<-
      as.list(1:(length(quant_splits)-1))
    
    for(chunck_num in seq_along(pw_comp_list)){
      start<-
        quant_splits[chunck_num]+1
      end<-
        quant_splits[chunck_num+1]
      pw_comp_list[[chunck_num]]<-
        start:end %>% 
        pw_comp_ids[.]
    }
    ############################################# 
    #pathway instances with name and list of components
    pw_inst<-
      pw_id %>% 
      selectInstances(biopax
                      ,id=.
                      ,includeReferencedInstances=FALSE)
    
    #get all referenced instances and convert them to biopax
    biopax_dt_list<-
      pw_comp_list %>% 
      lapply(selectInstances
             ,biopax = biopax
             ,includeReferencedInstances=TRUE) %>% 
      #append pw instances to each table
      lapply(rbind
             ,pw_inst) %>% 
      #convert to biopax
      lapply(biopax_from_dt) %>% 
      #remove unneeded components
      lapply(remove_deadend_refs)
    
    return(biopax_dt_list)
    
  }


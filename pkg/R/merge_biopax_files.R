merge_biopax_files<-
    function(source_dir=NULL
             ,change_ids=TRUE
             ,filename=NULL){
        
        #' @title
        #' Merge Multiple BioPAX Files into One File
        #' @description 
        #' Merge multiple BioPAX files into one BioPAX file, while preserving all unique relationships between pathway components.
        #' @param source_dir Directory with BioPAX files.
        #' @param change_ids Logical. Change ids to unify them throughout the files?
        #' @param filename Output filename (optional). If NULL (default), does not write the output to file.
        
        #' @author 
        #' Ivan Grishagin

        #in case a biopax source is in multiple files,
        #their tables are loaded and merged into one file
        require(rBiopaxParser
                ,quietly=TRUE)
        
        if (is.null(source_dir)){
            source_dir<-getwd()
        }
        
        #find owl directory and filenames in it
        owl.filename<-
            list.files(source_dir
                       ,pattern = "\\.owl"
                       ,full.names = TRUE)
        if(length(owl.filename)<1){
            stop("Biopax file has not been found in the following folder:\n"
                 ,owl.dirs[owl_num])
        }
     
        #read in the owl file(s), and merge their tables into one
        new_biopax_dt<-
            lapply(owl.filename
                   ,FUN=function(filename){
                       temp_biopax<-readBiopax(filename)
                       dt<-temp_biopax$dt
                       file_id<-
                           #get the filename without path
                           strsplit(filename
                                    ,split="\\/") %>% 
                           unlist %>%
                           .[length(.)] %>%
                           #if wiki - split by wp identifier
                           strsplit(split="\\_WP") %>% 
                           unlist %>%
                           .[length(.)] %>%
                           #remove owl extension
                           gsub("\\.owl","",.)
                       
                       #add id info to all pertaining entries
                       #get vector of ids
                       if(change_ids){
                           ids<-dt$id
                           clean_ref_ids<-
                               striphash(dt$property_attr_value)
                           
                           dt$id<-
                               paste(dt$id
                                     ,file_id
                                     #,sep="_pwref_")
                                     ,sep="_")
                           
                           dt$property_attr_value[clean_ref_ids %chin% ids]<-
                               dt$property_attr_value[clean_ref_ids %chin% ids] %>%
                               paste(file_id
                                     #,sep="_pwref_")
                                     ,sep="_")
                           
                       }
                       #change the dt in the temporary biopax object
                       temp_biopax$dt<-
                           dt
                       
                       #fix missing components
                       #get all non-referenced ids except for pathway one
                       #pathway ids
                       pwid<-
                           load_biopax_pathways(temp_biopax)$biopax.Pathway.ID
                       if(length(pwid)>1){
                           message("There's more than one pwid in the file "
                                   ,file_id
                                   ,"! Taking the first one with lowest number.")
                           pwid<-
                               sort(pwid)[1]
                       }
                       nonref_ids<-
                           temp_biopax$dt %>%
                           filter(!id %in% property_attr_value &
                                      !id %in% pwid) %>%
                           .$id %>%
                           unique
                       
                       if(length(nonref_ids)>0){
                           #add previously un-referenced ids as new components
                           temp_biopax<-
                               addPathwayComponents(temp_biopax
                                                    ,pwid
                                                    ,PATHWAY_COMPONENTS=nonref_ids)
                           message("Added pathway components "
                                   ,paste(nonref_ids
                                          ,collapse=", ")
                                   ," for the pathway "
                                   ,pwid
                                   ,".")
                       }
                      
                       return(temp_biopax$dt)
                   }) %>%
            do.call(rbind,.)
        #remove hashes if any
        new_biopax_dt$property_attr_value<-
            striphash(new_biopax_dt$property_attr_value)
            
        #create new biopax
        new_biopax<-
            new_biopax_dt %>% 
            biopax_from_dt(filename=filename)

        return(new_biopax)
    }
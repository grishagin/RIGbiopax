
############################### merge.biopax.files ############################### 
merge.biopax.files<-
    function(source_dir=NULL
             ,change_ids=TRUE
             ,write_to_file=FALSE){

        #in case a biopax source is in multiple files,
        #their tables are loaded and merged into one file
        require(rBiopaxParser
                ,quietly=TRUE)
        
        if (is.null(source_dir)){
            source_dir<-getwd()
        }
        
        #find owl directory and filenames in it
        owl.filename<-list.files(source_dir
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
                               rBiopaxParser:::striphash(dt$property_attr_value)
                           
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
                           load.biopax.pathways(temp_biopax)
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
                       
                       temp_biopax<-
                           addPathwayComponents(temp_biopax
                                                ,pwid
                                                ,PATHWAY_COMPONENTS=nonref_ids)
                       message("Added pathway components for the pathway "
                               ,pwid
                               ,".")
                      
                       return(dt)
                   }) %>%
            do.call(rbind,.)
        #create new biopax
        new_biopax<-
            createBiopax(level=3)
        new_biopax<-
            addBiopaxInstances(biopax=new_biopax
                               ,newInstancesDF=new_biopax_dt)

        
        if(write_to_file){
            biopax_filename<-
                paste(Sys.Date()
                      ,"combined_biopax.owl"
                      ,sep="_")
            writeBiopax_Rancho(biopax = new_biopax
                               ,file = biopax_filename
                               ,overwrite = TRUE
                               ,level=3)
        }

        return(new_biopax)
    }
############################### merge.biopax.files ############################### 
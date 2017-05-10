load_biopax_pathways<-
    function(owl_biopax=NULL){
        #' @title 
        #' Load BioPAX Pathways From BioPAX Object
        #' @description 
        #' Extracts all pathways from a BioPAX Object.
        #' @param owl_biopax A BioPAX object.
        #' 
        #' @author 
        #' Ivan Grishagin
        
        
        if(is.null(owl_biopax)){
            stop("Biopax object was not provided")
        }
        name_properties<-
            grep("name", owl_biopax$dt$property,ignore.case = TRUE, value = TRUE) %>% 
            unique() %>%
            tolower
        
        #pathway names are indicated in different biopax sources in a different manner
        #sometimes they are presented as standard name and name
        #and when they don't coincide, such "standardnames" are just GO processes
        if ("displayname" %in% name_properties) {
            #get instances of displayName of class pathway
            owl_biopax_instances_stname<-
                selectInstances(biopax=owl_biopax,
                                class="pathway",
                                property="displayName") 
            
            #subset a standardname df, select id and property columns
            pw_biopax<-
                owl_biopax_instances_stname %>%
                dplyr::select(id, name=property_value)
        } else if("standardname" %in% name_properties & 
                  "name" %in% name_properties ){
            #get instances of name and standardname of class pathway
            owl_biopax_instances_stname<-
                selectInstances(biopax=owl_biopax,
                                class="pathway",
                                property="standardName") 
            
            owl_biopax_instances_name<-
                selectInstances(biopax=owl_biopax,
                                class="pathway",
                                property="name")
            
            #instances of standrdname with same id as name
            matching_rows<-
                owl_biopax_instances_stname$id %in% owl_biopax_instances_name$id
            #subset a standardname df, select id and property columns
            pw_biopax<-
                owl_biopax_instances_stname[matching_rows,] %>%
                dplyr::select(id, name=property_value)
        } else if ("name" %in% name_properties) {
            #get instances of name of class pathway
            owl_biopax_instances_stname<-
                selectInstances(biopax=owl_biopax,
                                class="pathway",
                                property="name") 
            
            #subset a standardname df, select id and property columns
            pw_biopax<-
                owl_biopax_instances_stname %>%
                dplyr::select(id, name=property_value)
        } else {
            stop("load_biopax_pathways: check your name properties!")
        }
        
        #if pw_biopax is still empty, try default rBiopaxParser function
        if (nrow(pw_biopax)<1){
            pw_biopax<-
                listPathways(owl_biopax)
        }
        
        #convert names to lower case
        pw_biopax$name<-tolower(pw_biopax$name)
        
        #change colnum names
        colnames(pw_biopax)<-c("biopax.Pathway.ID","biopax.Pathway.Name")
        
        return(pw_biopax)
    }

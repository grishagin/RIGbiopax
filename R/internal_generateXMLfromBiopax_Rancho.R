internal_generateXMLfromBiopax_Rancho<-
    function (biopax
              ,namespaces = namespaces
              ,biopaxlevel=3) 
    {
        
        if (biopaxlevel==2){
            level_url<-"http://www.biopax.org/release/biopax-level2.owl"
        } else if (biopaxlevel==3){
            level_url<-"http://www.biopax.org/release/biopax-level3.owl"
        }
        
        #assemble the xml file header
        output_header<-
            '<?xml version="1.0" encoding="ISO-8859-1"?>\n<rdf:RDF '
        #add namespaces' definitions
        output_header<-
            paste0("xmlns:"
                   ,names(namespaces)
                   ,'="'
                   ,namespaces
                   ,'"') %>%
            paste(collapse="\n") %>%
            paste0(output_header
                   ,"\n"
                   ,.
                   ,">\n")
        #add ontology definitions
        output_header<-
            paste0('<owl:imports rdf:resource="'
                   ,level_url
                   ,'"/>\n'
                   ,'<rdfs:comment rdf:datatype="http://www.w3.org/2001/XMLSchema#string">'
                   ,'BioPAX output created '
                   ,date()
                   ,' using the Rancho-modified rBiopaxParser package.</rdfs:comment>\n'
            ) %>%
            paste0(output_header
                   ,'<owl:Ontology rdf:about="">\n'
                   ,.
                   ,'</owl:Ontology>')
        
        #add combined class/id id to the dataframe
        biopax$dt<-
            biopax$dt[,.(classid=paste0(class,id))] %>% 
            cbind(biopax$dt) %>% 
            .[order(classid)]
        
        #message("Combined ids prepared.")
        
        #open/close tag dataframe for each combined id
        instance_class_id<- 
            biopax$dt %>% 
            .[,.(classid
                 ,opentag=
                     paste0('<bp:'
                            ,class
                            ,' rdf:about="'
                            ,id
                            ,'">\n'
                     )
                 ,closetag=
                     paste0('</bp:'
                            ,class
                            ,'>'
                     ))] %>% 
            unique
        
        #now prepare tags for each instance in the biopax df
        #first, for those which have property value
        instance_df_pv<-
            biopax$dt %>%
            .[property_value!=""
              ,.(classid
                 ,tag=
                     paste0('<bp:'
                            ,property
                            ,' '
                            ,property_attr
                            ,'="'
                            ,property_attr_value
                            ,'">'
                            ,property_value
                            ,'</bp:'
                            ,property
                            ,'>\n'
                     ))]
            
        #... and then for those who don't
        instance_df_no_pv<-
            biopax$dt %>%
            .[property_value==""
              ,.(classid
                 ,tag=
                     paste0('<bp:'
                            ,property
                            ,' '
                            ,property_attr
                            ,'="'
                            ,property_attr_value
                            ,'"/>\n'))]
        
        #combine and turn in to a data table
        instance_df<-
            instance_df_pv %>%
            rbind(instance_df_no_pv) %>% 
            .[order(classid)]
        
        #set datatable key column
        setkey(instance_df,classid)
        
        
        #for each unique class/id combo find correponding referenced instances
        #and merge their corresponding tags into one line
        df_by_classid_kids<-
            instance_df[,.(kidstags=paste(tag,collapse="")),by=classid]
        
        #merge open tags with referenced instances tag strings and with close tags
        df_by_classid_full<-
            data.table(opentag=instance_class_id$opentag
                       ,kidstag=df_by_classid_kids$kidstags
                       ,closetag=instance_class_id$closetag) 
           
        #merge everything together into one string
        output<-
            df_by_classid_full %>% 
            .[,.(paste0(opentag
                        ,kidstag
                        ,closetag))] %>% 
            unlist %>% 
            c(output_header
              ,.
              ,"</rdf:RDF>")
    
        
        return(output)
    }

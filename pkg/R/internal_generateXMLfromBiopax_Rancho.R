internal_generateXMLfromBiopax_Rancho<-
    function (biopax
              ,filename=NULL
              ,namespaces = namespaces
              ,biopaxlevel=3
              ,verbose=TRUE) 
    {
        
        if (biopaxlevel==2){
            level_url<-"http://www.biopax.org/release/biopax-level2.owl"
        } else if (biopaxlevel==3){
            level_url<-"http://www.biopax.org/release/biopax-level3.owl"
        }
        
        #assemble the xml file header
        output<-
            '<?xml version="1.0"?>\n<rdf:RDF '
        #add namespaces' definitions
        output<-
            paste0("xmlns:"
                   ,names(namespaces)
                   ,'="'
                   ,namespaces
                   ,'"') %>%
            paste(collapse="\n") %>%
            paste0(output
                   ,"\n"
                   ,.
                   ,">\n")
        #add ontology definitions
        output<-
            paste0('<owl:imports rdf:resource="'
                   ,level_url
                   ,'"/>\n'
                   ,'<rdfs:comment rdf:datatype="http://www.w3.org/2001/XMLSchema#string">'
                   ,'BioPAX output created '
                   ,date()
                   ,' using the Rancho-modified rBiopaxParser package.</rdfs:comment>\n'
            ) %>%
            paste0(output
                   ,'<owl:Ontology rdf:about="">\n'
                   ,.
                   ,'</owl:Ontology>\n')
        
        #add combined class/id id to the dataframe
        biopax$dt<-
            biopax$dt %>%
            mutate(classid=paste0(class,id)) %>%
            arrange(classid)
        message("Combined ids prepared.")
        
        #open/close tag dataframe for each combined id
        instance_class_id<- 
            biopax$dt %>% 
            dplyr::select(class,id,classid) %>%
            unique %>%
            mutate(opentag=
                       paste0('<bp:'
                              ,class
                              ,' rdf:about="'
                              ,id
                              ,'">\n'
                       )
                   ,closetag=
                       paste0('</bp:'
                              ,class
                              ,'>\n'
                       )
            )
        
        #now prepare tags for each instance in the biopax df
        #first, for those with have property value
        instance_df_pv<-
            biopax$dt %>%
            filter(property_value!="") %>%
            mutate(tag=
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
                       )
            )
        #... and then for those who don't
        instance_df_no_pv<-
            biopax$dt %>%
            filter(property_value=="") %>%
            mutate(tag=
                       paste0('<bp:'
                              ,property
                              ,' '
                              ,property_attr
                              ,'="'
                              ,property_attr_value
                              ,'"/>\n'
                       )
            )
        #combine and turn in to a data table
        instance_df<-
            instance_df_pv %>%
            rbind.data.frame(instance_df_no_pv) %>%
            as.data.table
        
        message("All tags for interaction components and each instance were prepared.")
        
        #set datatable key column
        setkey(instance_df,classid)
        
        message("Adding instances' tags to referencing components' tags...")
        message("It should roughly "
                ,round(x = 0.002*nrow(instance_class_id)
                       ,digits = 0)
                ," seconds.")
        
        st<-Sys.time()
        #for each unique class/id combo find correponding referenced instances
        #and merge their corresponding tags into one line
        vect_by_classid_kids<-
            sapply(1:nrow(instance_class_id)
                   ,FUN=function(clid){
                       instance_df[instance_class_id$classid[clid],tag]
                   }) %>%
            lapply(paste,collapse="")
        et<-Sys.time()
        message("Tag addition is over. It took "
                ,round(x = difftime(et,st,units="secs")
                       ,digits = 0)
                ," seconds.")
        
        #merge open tags with referenced instances tag strings and with close tags
        vect_by_classid_full<-
            paste(instance_class_id$opentag
                  ,vect_by_classid_kids
                  ,instance_class_id$closetag
                  ,sep = "")
        #merge everything together into one string
        output<-
            vect_by_classid_full %>%
            paste(collapse="") %>%
            paste0(output
                   ,.
                   ,"</rdf:RDF>")
        
        if (is.null(filename)) {
            file<-paste0(Sys.Date()
                         ,"_generated_biopax.owl")
        }
        #write to file
        write(output
              ,file=filename)
        
        return(output)
    }

# internal_generateXMLfromBiopax_Rancho<-
#     function (biopax
#               ,namespaces = namespaces
#               ,verbose = TRUE
#               ,biopaxlevel=3) 
#     {
#         
#         if (biopaxlevel==2){
#             level_url<-"http://www.biopax.org/release/biopax-level2.owl"
#         } else if (biopaxlevel==3){
#             level_url<-"http://www.biopax.org/release/biopax-level3.owl"
#         }
#         
#         
#         d = suppressWarnings(XML::xmlTree("rdf:RDF"
#                                           ,namespaces = namespaces))
#         d$addNode("Ontology"
#                   ,namespace = "owl"
#                   ,attrs = c(`rdf:about` = "")
#                   ,close = FALSE)
#         d$addNode("imports"
#                   ,namespace = "owl"
#                   ,attrs = c(`rdf:resource` = level_url))
#         d$addNode("comment"
#                   ,paste("BioPAX output created"
#                          ,date()
#                          ,"using the Rancho-modified rBiopaxParser package.")
#                   ,namespace = "rdfs"
#                   ,attrs = c(`rdf:datatype` = "http://www.w3.org/2001/XMLSchema#string"))
#         
#         d$closeTag()
#         instanceList<- 
#             unique(biopax$dt[, list(class, id)])
#         count = 1
#         
#         #for (i in 1:dim(instanceList)[1]) {
#         for (i in 1:2000) {
#             instance = biopax$dt[class == instanceList[i]$class & 
#                                      id == instanceList[i]$id, ]
#             d$addNode(as.character(instance[1]$class)
#                       ,namespace = "bp"
#                       #,attrs = c(`rdf:ID` = as.character(instance[1]$id))
#                       #rdf:ID does not seem to work for SBML conversion -- yields empty document
#                       #in wiki, rdf:ID is replaced with rdf:about -- seems to work
#                       #inexplicably, original biocarta file (which has rdf:ID) works though
#                       ,attrs = c(`rdf:about` = as.character(instance[1]$id))
#                       ,close = FALSE)
#             
#             for (p in 1:dim(instance)[1]) {
#                 attrs = c(as.character(instance[p]$property_attr_value))
#                 names(attrs) = as.character(instance[p]$property_attr)
#                 
#                 if (nchar(as.character(instance[p]$property_value)) > 0) {
#                     d$addNode(as.character(instance[p]$property), 
#                               as.character(instance[p]$property_value)
#                               ,namespace = "bp"
#                               ,attrs = attrs)
#                 }
#                 else {
#                     d$addNode(as.character(instance[p]$property)
#                               ,namespace = "bp"
#                               ,attrs = attrs)
#                 }
#             }
#             d$closeTag()
#             if (verbose) {
#                 if (i%%1000 == 0) 
#                     message(paste0("INFO: Wrote instance nr "
#                                    ,i
#                                    ," of "
#                                    ,dim(instanceList)[1]
#                                    ," with id "
#                                    ,instance$id[1]
#                                    ,"."))
#                 count = count + 1
#             }
#         }
#         d
#     }

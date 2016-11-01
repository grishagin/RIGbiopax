internal_generateXMLfromBiopax_Rancho<-
    function (biopax
              ,namespaces = namespaces
              ,verbose = TRUE
              ,biopaxlevel=3) 
    {
        
        if (biopaxlevel==2){
            level_url<-"http://www.biopax.org/release/biopax-level2.owl"
        } else if (biopaxlevel==3){
            level_url<-"http://www.biopax.org/release/biopax-level3.owl"
        }
        
        
        d = suppressWarnings(XML::xmlTree("rdf:RDF"
                                          ,namespaces = namespaces))
        d$addNode("Ontology"
                  ,namespace = "owl"
                  ,attrs = c(`rdf:about` = "")
                  ,close = FALSE)
        d$addNode("imports"
                  ,namespace = "owl"
                  ,attrs = c(`rdf:resource` = level_url))
        d$addNode("comment"
                  ,paste("BioPAX output created"
                         ,date()
                         ,"using the Rancho-modified rBiopaxParser package.")
                  ,namespace = "rdfs"
                  ,attrs = c(`rdf:datatype` = "http://www.w3.org/2001/XMLSchema#string"))
        
        d$closeTag()
        instanceList<- 
            unique(biopax$dt[, list(class, id)])
        count = 1
        
        for (i in 1:dim(instanceList)[1]) {
            instance = biopax$dt[class == instanceList[i]$class & 
                                     id == instanceList[i]$id, ]
            d$addNode(as.character(instance[1]$class)
                      ,namespace = "bp"
                      #,attrs = c(`rdf:ID` = as.character(instance[1]$id))
                      #rdf:ID does not seem to work for SBML conversion -- yields empty document
                      #in wiki, rdf:ID is replaced with rdf:about -- seems to work
                      #inexplicably, original biocarta file (which has rdf:ID) works though
                      ,attrs = c(`rdf:about` = as.character(instance[1]$id))
                      ,close = FALSE)
            
            for (p in 1:dim(instance)[1]) {
                attrs = c(as.character(instance[p]$property_attr_value))
                names(attrs) = as.character(instance[p]$property_attr)
                
                if (nchar(as.character(instance[p]$property_value)) > 0) {
                    d$addNode(as.character(instance[p]$property), 
                              as.character(instance[p]$property_value)
                              ,namespace = "bp"
                              ,attrs = attrs)
                }
                else {
                    d$addNode(as.character(instance[p]$property)
                              ,namespace = "bp"
                              ,attrs = attrs)
                }
            }
            d$closeTag()
            if (verbose) {
                if (i%%1000 == 0) 
                    message(paste("INFO: Wrote instance nr"
                                  ,i
                                  ,"of"
                                  ,dim(instanceList)[1]
                                  ,"with id"
                                  ,instance$id[1]
                                  ,"."))
                count = count + 1
            }
        }
        d
    }

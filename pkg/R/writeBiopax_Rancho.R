writeBiopax_Rancho<-
        function (biopax
                  ,file = NULL
                  ,verbose = TRUE
                  ,overwrite = FALSE
                  ,biopaxlevel=3
                  ,namespaces=NULL) {
            
                if (!biopaxlevel %in% c(2,3)){
                    stop("writeBiopax_Rancho: Incorrect biopax level specified. Aborting.")
                }

            
                if (is.null(namespaces)){
                    if (biopaxlevel==2){
                        message("Using default namespaces for biopax level 2.")
                        namespaces<-list(rdf = "http://www.w3.org/1999/02/22-rdf-syntax-ns#", 
                                         bp = "http://www.biopax.org/release/biopax-biopaxlevel2.owl#", 
                                         rdfs = "http://www.w3.org/2000/01/rdf-schema#", owl = "http://www.w3.org/2002/07/owl#", 
                                         xsd = "http://www.w3.org/2001/XMLSchema#")
                    } else if (biopaxlevel==3){
                        message("Using default namespaces for biopax level 3.")
                        namespaces<-list(rdf = "http://www.w3.org/1999/02/22-rdf-syntax-ns#", 
                                         bp = "http://www.biopax.org/release/biopax-biopaxlevel3.owl#", 
                                         rdfs = "http://www.w3.org/2000/01/rdf-schema#", owl = "http://www.w3.org/2002/07/owl#", 
                                         xsd = "http://www.w3.org/2001/XMLSchema#")
                    }
                } 
                
                
                if (file.exists(file) & !overwrite) {
                        stop(paste("Error: File ", file, " already exists.", 
                                   sep = ""))
                }
                
                .checkValidity_Rancho(biopax)
                
                d <-
                    .internal_generateXMLfromBiopax_Rancho(biopax
                                                           ,namespaces
                                                           ,verbose = verbose
                                                           ,biopaxlevel = biopaxlevel)
                if (is.null(file)) {
                        file=paste0(Sys.Date()
                                    ,"_generated_biopax.owl")
                }
                XML::saveXML(d$value()
                             ,file = file
                             ,encoding = "ISO-8859-1")
        }
writeBiopax_Rancho<-
        function (biopax
                  ,file = NULL
                  ,verbose = TRUE
                  ,overwrite = FALSE
                  ,level=3
                  ,namespaces=NULL) {

                if (level==2){
                        namespaces<-list(rdf = "http://www.w3.org/1999/02/22-rdf-syntax-ns#", 
                                         bp = "http://www.biopax.org/release/biopax-level2.owl#", 
                                         rdfs = "http://www.w3.org/2000/01/rdf-schema#", owl = "http://www.w3.org/2002/07/owl#", 
                                         xsd = "http://www.w3.org/2001/XMLSchema#")
                } else if (level==3){
                        namespaces<-list(rdf = "http://www.w3.org/1999/02/22-rdf-syntax-ns#", 
                                         bp = "http://www.biopax.org/release/biopax-level3.owl#", 
                                         rdfs = "http://www.w3.org/2000/01/rdf-schema#", owl = "http://www.w3.org/2002/07/owl#", 
                                         xsd = "http://www.w3.org/2001/XMLSchema#")
                } else if (is.null(level)){
                        if(is.null(namespaces)){
                                stop("Since you have provided a biopax level NULL, you have to specify namespaces. Stopping.")
                        }
                }
                
                
                if (file.exists(file) & !overwrite) {
                        stop(paste("Error: File ", file, " already exists.", 
                                   sep = ""))
                }
                rBiopaxParser:::checkValidity(biopax)
                d = rBiopaxParser:::internal_generateXMLfromBiopax(biopax, namespaces, verbose = verbose)
                if (is.null(file)) {
                        file=paste0(Sys.Date()
                                    ,"_generated_biopax.owl")
                }
                XML::saveXML(d$value()
                             ,file = file)
        }
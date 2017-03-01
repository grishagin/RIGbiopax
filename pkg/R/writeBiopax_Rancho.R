writeBiopax_Rancho<-
    function (biopax
              ,filename = "default"
              ,overwrite = FALSE
              ,biopaxlevel=3
              ,namespaces=NULL) {
        
        if (!biopaxlevel %in% c(2,3)){
            stop("writeBiopax_Rancho: Incorrect biopax level specified. Aborting.")
        }
        
        
        if (is.null(namespaces)){
            if (biopaxlevel==2){
                message("Using default namespaces for biopax level 2.")
                namespaces<-
                    list(rdf = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
                         ,bp = "http://www.biopax.org/release/biopax-level2.owl#" 
                         ,rdfs = "http://www.w3.org/2000/01/rdf-schema#"
                         ,owl = "http://www.w3.org/2002/07/owl#"
                         ,xsd = "http://www.w3.org/2001/XMLSchema#")
            } else if (biopaxlevel==3){
                message("Using default namespaces for biopax level 3.")
                namespaces<-
                    list(rdf = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
                         ,bp = "http://www.biopax.org/release/biopax-level3.owl#"
                         ,rdfs = "http://www.w3.org/2000/01/rdf-schema#"
                         ,owl = "http://www.w3.org/2002/07/owl#"
                         ,xsd = "http://www.w3.org/2001/XMLSchema#")
            }
        } 
        
        
        if (file.exists(filename) & !overwrite) {
            stop(paste("Error: File ", filename, " already exists.", 
                       sep = ""))
        }
        
        RIGbiopax:::internal_checkValidity_Rancho(biopax)
        
        output <-
            internal_generateXMLfromBiopax_Rancho(biopax = biopax
                                                  ,filename=filename
                                                  ,namespaces = namespaces
                                                  ,biopaxlevel = biopaxlevel
                                                  ,verbose=verbose)
        if (filename=="default") {
            filename<-
                paste0(Sys.Date()
                       ,"_generated_biopax.owl")
        }
        #write to file
        writeLines(output
              ,con = filename
              ,useBytes = TRUE)
        return()
    }

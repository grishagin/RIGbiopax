writeBiopax_Rancho<-
  function (biopax
            ,filename = "default"
            ,overwrite = FALSE
            ,biopaxlevel=3
            ,namespaces=NULL
            ,verbose=TRUE) {
    #' @export
    #' @title
    #' Write a BioPAX Object into File
    #' @description 
    #' Converts a BioPAX object into an XML object, and writes it to file.
    #' @param biopax A BioPAX object.
    #' @param filename (string) Output file name.
    #' @param overwrite (boolean) Overwrite the file if already exists?
    #' @param namespaces (list) A list of names strings required 
    #' to construct an XML file. If \code{NULL} will default to values 
    #' appropriate for the corresponding BioPAX level.
    #' @param verbose (boolean) Show all messages?
    
    #' @author 
    #' Ivan Grishagin
    
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
    
    if (filename=="default") {
      filename<-
        paste0(Sys.Date()
               ,"_generated_biopax.owl")
    }
    
    if(verbose){
      message("Writing BioPAX object to file "
              ,filename
              ,".\nIt should take "
              ,round(x = 1.55e-5*nrow(biopax$dt)
                     ,digits = 0)
              ," seconds.")
    }
    
    st<-Sys.time()
    output<-
      internal_generateXMLfromBiopax_Rancho(biopax = biopax
                                            ,namespaces = namespaces
                                            ,biopaxlevel = biopaxlevel)
    #write to file
    writeLines(output
               ,con = filename
               ,useBytes = TRUE)
    
    et<-Sys.time()
    
    if(verbose){
      message("Writing BioPAX object is complete. It took "
              ,round(x = difftime(et,st,units="secs")
                     ,digits = 2)
              ," seconds.")
    }
    return()
  }

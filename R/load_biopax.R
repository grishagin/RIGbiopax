load_biopax<-
  function(source_name=NULL
           ,source_dir=NULL){
    
    #' @export
    #' @title 
    #' Load BioPAX from Environment or OWL or Workspace File
    #' @description 
    #' Tries to load BioPAX object in several ways: searches for a BioPAX object in the environment, 
    #' and prompts if more than one object or no object has been found. Otherwise, looks for an \code{*.owl} file in the specified directory. 
    #' Returns a BioPAX object.
    #' @param source_name Name of a BioPAX source.
    #' @param source_dir Directory with BioPAX files (optional).
    #' 
    #' @author 
    #' Ivan Grishagin
    
    #get all obj names in global environment
    objects<-
      ls(envir = .GlobalEnv)
    
    #if params for reading biopax have not been specified...
    if (is.null(source_name) | is.null(source_dir)){
      
      #check for an object amongst loaded objects
      is_biopax<-
        sapply(objects,function(obj) {
          "biopax" %in% class(get(obj))   
        }) 
      if(sum(is_biopax)<1){
        #if no biopax object found in the environment,
        #prompt to choose workspace or biopax file
        message("Biopax object was not found in the environment!")
        
        #... ask for workspace or biopax file
        filetoload<-
          invisible(tcltk::tk_choose.files(caption = "Choose *.owl or *.Rdata file to read/load."
                                           ,multi=FALSE
                                           ,filters=matrix(c("BioPAX"
                                                             ,"Workspace"
                                                             ,".owl"
                                                             ,".Rdata")
                                                           ,ncol=2)))
        #check if in either case the file can be loaded/read
        trytoloadws<-
          suppressWarnings(try(load(filetoload),silent=TRUE))
        trytoreadbp<-
          suppressWarnings(try(biopax<-
                                 readBiopax(filetoload),silent=TRUE))
        
        #if not, inform
        if(trytoloadws!=TRUE & trytoreadbp!=TRUE){
          message("Could not interpret the file as either a Workspace or a BioPAX!")
        }
        
        #successful or not, try loading biopax again
        biopax<-
          load_biopax()
        
      } else if(sum(is_biopax)>1){
        #if more than one biopax, 
        #open a radiobutton picker form
        message("More than one biopax object present in environment!")
        message("Pick the desired object name in a popup window")
        
        biopax<-
          get(tkradio_from_vect(vect = subset(objects
                                              ,is_biopax)
                                ,window_title = "Pick the desired biopax object"))
        
      } else {
        #if just one object found,
        #return it as biopax
        biopax<-
          get(subset(objects
                     ,is_biopax))
      }
      
    } else {
      #if both file name and dir name have been provided
      #try to load the file
      biopax<-
        load_biopax_file(source_name=source_name
                         ,source_dir=source_dir)
    }
    
    return(biopax)
  }

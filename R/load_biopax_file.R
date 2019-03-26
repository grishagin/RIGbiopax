load_biopax_file <-
    function(source_name=NULL
             ,source_dir=NULL){
        #' @title 
        #' Load BioPAX from OWL File
        #' @description 
        #' Tries to load BioPAX object from an OWL file.
        #' @param source_name Name of a BioPAX source.
        #' @param source_dir Directory with BioPAX files (optional).   
        #'  
        #' @author 
        #' Ivan Grishagin
        
        
        require(rBiopaxParser)
        if(is.null(source_name)){
            stop("Provide source_name!")
        }
        if (is.null(source_dir)){
            source_dir<-getwd()
        }
        
        #find owl directory
        owl_dirs<-list.dirs(source_dir
                            ,recursive = FALSE)
        owl_num<-agrep(source_name
                       ,owl_dirs
                       ,max.distance = 0.2)
        
        owl.filename<-list.files(owl_dirs[owl_num]
                                 ,pattern = "\\.owl"
                                 ,full.names = TRUE)
        #check if only one filename has been fetched
        if(length(owl.filename)!=1){
            stop("There's not just one biopax file in the following folder:\n",owl_dirs[owl_num])
        }
        #read in the owl file
        owl_biopax<-readBiopax(owl.filename)
    }

load.biopax <-
function(source_name=NULL
             ,source_dir=NULL){
        #get all objects from parent environment
        objects<-
            ls(envir = parent.env(env = environment()))
        #any biopax files in env?
        is_present_biopax<-
            sapply(objects,function(obj) {
                "biopax" %in% class(get(obj))   
            }) 
        
        if(sum(is_present_biopax)<1){
            biopax<-
                load.biopax.file(source_name=source_name
                                 ,source_dir=source_dir)
        } else if(sum(is_present_biopax)>1){
            stop("More than one biopax object present in environment already")
        } else {
            biopax<-
                get(subset(objects,is_present_biopax))
        }
        return(biopax)
    }

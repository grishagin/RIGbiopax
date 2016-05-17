load.biopax.file <-
function(source_name=NULL
             ,source_dir=NULL){
        require(rBiopaxParser)
        if(is.null(source_name)){
            stop("Provide source_name!")
        }
        if (is.null(source_dir)){
            source_dir<-getwd()
        }
        
        #find owl directory
        owl.dirs<-list.dirs(source_dir
                            ,recursive = FALSE)
        owl_num<-agrep(source_name
                       ,owl.dirs
                       ,max.distance = 0.2)
        
        owl.filename<-list.files(owl.dirs[owl_num]
                                 ,pattern = "\\.owl"
                                 ,full.names = TRUE)
        #check if only one filename has been fetched
        if(length(owl.filename)!=1){
            stop("There's not just one biopax file in the following folder:\n",owl.dirs[owl_num])
        }
        #read in the owl file
        owl_biopax<-readBiopax(owl.filename)
    }

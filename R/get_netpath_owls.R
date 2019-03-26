get_netpath_owls<-
    function(dir=getwd()){
        #' @title
        #' Download NetPath BioPAX Files
        #' @description 
        #' Downloads NetPath BioPAX files in the *.owl format.
        #' @param dir Directory, where the files will be saved to in a new subdirectory \code{NetPath_owls}.
        #' @author 
        #' Ivan Grishagin
        #' 
        #' 
        
        require(XML)
        require(RCurl)
        
        #parse html contents
        netpath.HTML<-htmlParse(getURL("http://www.netpath.org/download/BioPAX"))
        
        #get a set of nodes which are links in a table
        hrefs.list<-getNodeSet(netpath.HTML,"//table/..//a/@href") %>%
            as.character
        
        #get a list of owl urls
        owl.urls.list<-grepl(".owl"
                             ,hrefs.list
                             ,ignore.case=TRUE) %>% 
            hrefs.list[.]
        
        #download files
        dir.create(file.path(dir
                             ,"NetPath_owls")
                   ,showWarnings = FALSE)
        
        for (index in 1:length(owl.urls.list)){
            file_name<-
                strsplit(owl.urls.list[[index]]
                         ,split="/") %>%
                unlist %>%
                .[length(.)] %>%
                file.path(dir
                          ,"NetPath_owls"
                          ,.)
            
            download.file(url=owl.urls.list[[index]]
                          ,destfile=file_name
                          ,quiet = TRUE)
        }
        
    }
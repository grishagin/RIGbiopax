clean_biopax_property_value<-
    function(biopax){
        #function cleans up property value column
        #1) removes html tags
        #2) removes eol chars
        #3) replaces entire uniprot (?) entries that start with "Recname:" 
        #with the recommended names
        
        biopax$dt$property_value<-
            biopax$dt$property_value %>%
            gsub(pattern="\n|<[^=]*?>|<[^-]*?>|<[^[:space:]]*?>"
                 ,replacement = ""
                 ,.)
        
        biopax$dt$property_value<-
            biopax$dt$property_value %>%
            gsub(pattern="^RecName: Full=(.*?);.*$"
                 ,replacement = "\\1"
                 ,.)
        return(biopax)
    }
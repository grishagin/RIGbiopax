clean_biopax<-
    function(biopax){
        
        #' @title
        #' Clean BioPAX Object
        #' @description 
        #' 1) Strips hashes from property_value column.
        #' 2) Cleans a BioPAX object from html tags and garbled symbols. 
        #' 3) Unifies representation of some classes, common properties and property values.
        #' 4) Removes deadend references.
        #' @details 
        #'
        #' @param biopax BioPAX object.
       
        #' @author 
        #' Ivan Grishagin
        
        ################################ strip hashes from property value column
        biopax$dt$property_attr_value<-
            biopax$dt$property_attr_value %>% 
            striphash
        
        ################################ clean up property value column
        #1) removes html tags
        #2) removes eol chars
        #3) replaces entire uniprot (?) entries that start with "Recname:" with the recommended names
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
        biopax$dt$property_value[biopax$dt$property_value=="NA"]<-
            "UNKNOWN"
        
        ################################ garbled symbols
        dict<-
            data.table(from=c("\u201A\u00e0\u00f6\u221a\u00a9\u00ac\u00a8\u00ac\u00b5"
                              ,"\u00c3\u008e\u00c2\u00b2"
                              ,"\u00c3\u009f","&gt;","&apos;","&"
                              ,"\u00e2\u0080\u009c"
                              ,"\u00e2\u0080?"
                              ,"\u00c3\u00a2\u00e2\u0082\u00ac\u00e2\u0080\u009c"
                              ,"\u00c3\u00a2\u00e2\u0082\u00ac\u00e2"
                              ,"\u00e2\u0093"
                              ,"\u00e2\u0094"
                              ,"\u00c3\u00a2\u00e2\u0082\u00ac\u00e2\u0084\u00a2"
                              ,"\u00c3\u00a2\u00e2\u0082\u00ac\u00c2\u00b2"
                              ,"\u00c3\u00a2\u00e2\u0082\u00ac\u00c5\u0093"
                              ,"\u00c3\u00a2\u00e2\u0082\u00ac\u00c2?"
                              ,"\u0093"
                              ,"\u0094")
                       ,to=c("\u03b5"
                             ,"\u03b2"
                             ,"\u03b2",">","'","and"
                             ,""
                             ,""
                             ,"-"
                             ,"-"
                             ,"-"
                             ,"-"
                             ,"'"
                             ,"'"
                             ,"\""
                             ,"\""
                             ,"\""
                             ,"\""))
        
        #clean biopax object from utf tags
        #find affected rows
        badsymbol_rows<-
            biopax$dt[property_value!=""]$property_value %>% 
            grepl("[^[:graph:][:space:]]"
                  ,.) 
        
        biopax$dt[property_value!=""][badsymbol_rows]$property_value<-
            biopax$dt[property_value!=""][badsymbol_rows]$property_value %>%
            mgsub(pattern = dict$from
                  ,replacement = dict$to
                  ,text.var = .)
        
        ################################ class
        class_dict<-
            data.table(from=c("Rna"
                              ,"RnaReference"
                              ,"Dna"
                              ,"DnaReference")
                       ,to=c("RNA"
                             ,"RNAReference"
                             ,"DNA"
                             ,"DNAReference"))
                              
        biopax$dt$class<-
            biopax$dt$class %>% 
            mapvalues(from=class_dict$from
                      ,to=class_dict$to
                      ,warn_missing = FALSE)
        
        ################################ property
        property_dict<-
            data.table(from=c("dispName")
                       ,to=c("displayName"))
                       
        biopax$dt$property<-
            biopax$dt$property %>% 
            mapvalues(from=property_dict$from
                      ,to=property_dict$to
                      ,warn_missing = FALSE)
        
        ################################ property_value
        property_value_dict<-
            data.table(from=c("LEFT_to_RIGHT")
                       ,to=c("LEFT-to-RIGHT"))
        
        biopax$dt$property_value<-
            biopax$dt$property_value %>% 
            mapvalues(from=property_value_dict$from
                      ,to=property_value_dict$to
                      ,warn_missing = FALSE)
        
        #unify database names
        db_dict<-
            system.file("extdata"
                        ,"unify_db_names.xlsx"
                        ,package="RIGbiopax") %>% 
            read_excel_astext(sheet = 1)
        
        biopax$dt$property_value<-
            biopax$dt$property_value %>% 
            mapvalues(from=db_dict$original_db_name
                      ,to=db_dict$new_db_name
                      ,warn_missing = FALSE)
        
        ################################ remove deadend references
        biopax<-
            biopax %>% 
            remove_deadend_refs
        
        
        return(biopax)
    }
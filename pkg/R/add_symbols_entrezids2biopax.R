add_symbols_entrezids2biopax<-
    function(biopax){
        
        #' @title 
        #' Annotate a BioPAX Object with Gene ID and Gene Symbol
        #' @description 
        #' Adds Gene IDs and Symbols to all physical entities, where possible.
        #' @param biopax A BioPAX object.
        #' 
        #' @author 
        #' Ivan Grishagin
        
        #find all current dbid pairs 
        dbid_df<-
            biopax$dt[property %in% c("db"
                                      ,"id")
                      ,.(dbid_db=property_value[property=="db"]
                         ,dbid_id=property_value[property=="id"])
                      ,by=id] %>% 
            unique
        
        #find xrefs referencing these ids
        xref_df<-
            biopax$dt[property_attr_value %in% dbid_df$id
                      ,.(xref_id=id
                         ,pav=property_attr_value)] %>% 
            unique
        
        #match them up and add to main df
        dbid_df<-
            dbid_df[match(xref_df$pav
                      ,id)] %>% 
            cbind(xref_df[,.(xref_id)]
                  ,.)
        
        #find physent referencing these xrefs
        physent_df<-
            biopax$dt[property_attr_value %in% dbid_df$xref_id
                      ,.(physent_id=id
                         ,physent_class=class
                         ,xref_pav=property_attr_value)] %>% 
            unique

        #match them up and add to main df
        dbid_df<-
            dbid_df[match(physent_df$xref_pav
                          ,xref_id)] %>% 
            cbind(physent_df[,.(physent_id
                                ,physent_class)]
                  ,.)
        
        #annotate the data table with entrez ids and symbols
        dbid_df<-
            dbid_df %>% 
            add_symbols_entrezids_mult_keytypes(col_names=
                                                    c(id_col = "dbid_id"
                                                      ,type_col = "dbid_db"
                                                      ,entrez_col = "entrezgene"
                                                      ,symbol_col = "symbol"))
        
        #prepare all actual data tables that connect
        #physical entities to the entrez ids and symbols
        dTable_list<-list()
        dTable_list$entref<-
            dbid_df[,.(class=physent_class
                       ,id=physent_id
                       ,property="entityReference"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=
                           paste0(xref_id
                                  ,"entrez")
                       ,property_value="")]
        
        dTable_list$symbol<-
            dbid_df[,.(class=
                           paste0(physent_class
                                  ,"Reference")
                       ,id=
                           paste0(xref_id
                                  ,"entrez")
                       ,property="name"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=symbol)]
        
        if(nrow(dbid_df[dbid_id!=entrezgene])>0){
            dTable_list$xref<-
                dbid_df[,.(class=
                               paste0(physent_class[dbid_id!=entrezgene]
                                      ,"Reference")
                           ,id=
                               paste0(xref_id[dbid_id!=entrezgene]
                                      ,"entrez")
                           ,property="xref"
                           ,property_attr="rdf:resource"
                           ,property_attr_value=
                               paste0(id[dbid_id!=entrezgene]
                                      ,"entrez")
                           ,property_value="")]
            dTable_list$db<-
                #do not take those values, where dbid_id==entrezid
                dbid_df[,.(class="RelationshipXref"
                           ,id=
                               paste0(id[dbid_id!=entrezgene]
                                      ,"entrez")
                           ,property="db"
                           ,property_attr="rdf:datatype"
                           ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                           ,property_value="entrezgene")]
            
            dTable_list$id<-
                #do not take those values, where dbid_id==entrezid
                dbid_df[,.(class="RelationshipXref"
                           ,id=
                               paste0(id[dbid_id!=entrezgene]
                                      ,"entrez")
                           ,property="id"
                           ,property_attr="rdf:datatype"
                           ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                           ,property_value=entrezgene[dbid_id!=entrezgene])]
        }
        
        
        #merge them and add them to biopax$dt
        biopax$dt<-
            dTable_list %>% 
            do.call(rbind
                    ,.) %>% 
            unique %>% 
            rbind(biopax$dt)

        
        return(biopax)
    }
gene_bag_to_biopax<-
    function(gene_df=NULL
             ,filename=NULL
             ,cols=list(pwid="Pathway.ID"
                        ,pwname="Pathway.Name"
                        ,genedisplayname="Gene.Symbol"
                        ,genename="Gene.Name"
                        ,genesymbol="Gene.Symbol"
                        ,geneid="Gene.ID")
             ){
        #ensure that gene_df is a dataframe
        gene_df<-
            gene_df %>%
            as.data.frame
        #declare list to store all biopax-style data table
        dTable_list<-
            list()
        #add unique gene name/display name combinations to the dataframe
        gene_df$cmbn<-
            paste0(gene_df[,cols$genedisplayname]
                   ,gene_df[,cols$genename])
        
        #prepare component ids
        comp_ids<-
            paste(gene_df[,cols$pwid]
                  ,unlist(lapply(unique(gene_df[,cols$pwid])
                                 ,FUN = function(id){
                                     #rows for given pw id
                                     rows<-
                                         (gene_df[,cols$pwid] %in% id)
                                     #prepare indices for unique gene name/displayname
                                     #combinations for given pw
                                     ind<-
                                         match(gene_df$cmbn[rows]
                                               ,unique(gene_df$cmbn[rows]))
                                     return(ind)
                                 }))
                  ,"c"
                  ,sep = "_")
        
        
        
        #references to the component names/xrefs/etc
        xref_ids<-
            paste0(comp_ids
                   ,"x"
                   ,internal_seq_along_find_reps(comp_ids))
        #references to db and id names
        db_ids<-
            paste0(xref_ids
                   ,"x")
        
        #declare a biopax-style data table and fill pathway components
        dTable_list$dt_comp_ids<-
            data.frame(class="Pathway"
                       ,id=gene_df[,cols$pwid]
                       ,property="pathwayComponent"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=comp_ids
                       ,property_value=""
            )
        #pathway names
        dTable_list$dt_pw_names<-
            data.frame(class="Pathway"
                       ,id=gene_df[,cols$pwid]
                       ,property="name"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df[,cols$pwname]
            ) %>%
            unique
        #protein display names
        dTable_list$dt_prot_displaynames<-
            data.frame(class="Protein"
                       ,id=comp_ids
                       ,property="displayName"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df[,cols$genedisplayname]
            )
        #protein names
        dTable_list$dt_prot_names<-
            data.frame(class="Protein"
                       ,id=comp_ids
                       ,property="name"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df[,cols$genename]
            )
        
        #protein entity references
        dTable_list$dt_prot_refs<-
            data.frame(class="Protein"
                       ,id=comp_ids
                       ,property="entityReference"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=xref_ids
                       ,property_value=""
            )
        
        #protein gene symbol
        dTable_list$dt_prot_refs_sym<-
            data.frame(class="ProteinReference"
                       ,id=xref_ids
                       ,property="name"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df[,cols$genesymbol]
            )
        
        #protein entity references to the gene 
        dTable_list$dt_prot_refs_x<-
            data.frame(class="ProteinReference"
                       ,id=xref_ids
                       ,property="xref"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=db_ids
                       ,property_value=""
            )
        
        #protein id db
        dTable_list$dt_prot_refs_x_db<-
            data.frame(class="RelationshipXref"
                       ,id=db_ids
                       ,property="db"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value="entrezgene"
            )
        
        #protein gene id
        dTable_list$dt_prot_refs_x_id<-
            data.frame(class="RelationshipXref"
                       ,id=db_ids
                       ,property="id"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df[,cols$geneid]
            )
        dTable<-
            do.call(rbind
                    ,dTable_list) 
        
        biopax<-
            biopax_from_dt(dTable
                           ,filename=filename)
        
    }



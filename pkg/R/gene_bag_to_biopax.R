gene_bag_to_biopax<-
    function(gene_df=NULL
             ,filename=NULL){
        mandatory_cols<-
            c("Pathway.ID"
              ,"Pathway.Name"
              ,"Gene.Symbol"
              ,"Gene.ID")
        if (is.null(colnames(gene_df))){
            stop("gene_bag_to_biopax: gene_df does not appear to have columns!")
        } else if (!all(colnames(gene_df) %in% mandatory_cols)){
            stop("gene_bag_to_biopax: gene_df does not have all mandatory columns!")
        }
        #declare list to store all biopax-style data table
        dTable_list<-
            list()
        #prepare component ids
        comp_ids<-
            paste0(gene_df$Pathway.ID
                   ,unlist(lapply(unique(gene_df$Pathway.ID)
                                  ,FUN = function(id){
                                      1:sum(gene_df$Pathway.ID %in% id)
                                  }))
                   ,"_c")
        xref_ids<-
            paste0(gene_df$Pathway.ID
                   ,"_cx")
        db_ids<-
            paste0(gene_df$Pathway.ID
                   ,"_cxx")
        #declare a biopax-style data table and fill pathway components
        dTable_list$dt_comp_ids<-
            data.frame(class="Pathway"
                       ,id=gene_df$Pathway.ID
                       ,property="pathwayComponent"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=comp_ids
                       ,property_value=""
            )
        #pathway names
        dTable_list$dt_pw_names<-
            data.frame(class="Pathway"
                       ,id=gene_df$Pathway.ID
                       ,property="name"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df$Pathway.Name
            ) %>%
            unique
        #protein names
        dTable_list$dt_prot_names<-
            data.frame(class="Protein"
                       ,id=comp_ids
                       ,property="name"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df$Gene.Symbol
            )
        
        #protein entity references
        dTable_list$dt_prot_refs<-
            data.frame(class="ProteinReference"
                       ,id=comp_ids
                       ,property="entityReference"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=xref_ids
                       ,property_value=""
            )
        #more protein entity references 
        dTable_list$dt_prot_refs_x<-
            data.frame(class="RelationshipXref"
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
                       ,property_value=gene_df$Gene.ID
            )
        dTable<-
            do.call(rbind
                    ,dTable_list)
        
        biopax<-
            biopax_from_dt(dTable
                           ,filename=filename)
        
    }



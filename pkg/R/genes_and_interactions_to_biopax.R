genes_and_interactions_to_biopax<-
    function(gene_df=NULL
             ,filename=NULL
             ,cols=list(pwid=
                            "toxdb.Pathway.ID"
                        ,srcpwname=
                            "source.Pathway.Name"
                        
                        ,first_dispname=
                            "first_prot_sym"
                        ,first_name=
                            "first_name"
                        ,first_genesym=
                            "first_gene_sym"
                        ,first_geneid=
                            "first_gene_id"
                        
                        ,second_dispname=
                            "second_prot_sym"	
                        ,second_name=
                            "second_name"
                        ,second_genesym=
                            "second_gene_sym"
                        ,second_geneid=
                            "second_gene_id"
                        
                        ,pmid_author=
                            "pmid_author"
                        
                        ,rxn_class=
                            "reaction_class"
                        ,ctrl_class=
                            "control_class"
                        ,ctrl_type=
                            "control_type"
                        ,second_left_right=
                            "second_left_right"
                        )
    ){
        #ensure that gene_df is a dataframe
        gene_df_orig<-
            gene_df %>%
            as.data.frame
        #replace all text "NA" with real NA values
        gene_df_orig[gene_df_orig %in% "NA"]<-
            NA

        #first component left or right property
        #NA by default, and left only if 
        #1. the second component is left or right and
        #2. it's not a control interaction
        gene_df_orig$first_left_right<-
            NA
        gene_df_orig$first_left_right[!is.na(gene_df_orig[,cols$second_left_right]) &
                                          is.na(gene_df_orig[,cols$ctrl_class])]<-
            "left"
        
        gene_df_orig$first_ctrl_prop<-
            NA
        gene_df_orig$first_ctrl_prop[!is.na(gene_df_orig[,cols$ctrl_class])]<-
            "controller"
        gene_df_orig$second_ctrl_prop<-
            NA
        gene_df_orig$second_ctrl_prop[!is.na(gene_df_orig[,cols$ctrl_class])]<-
            "controlled"
        
        #make up control component ids
        gene_df_orig$ctrl_id<-
            paste(gene_df_orig[,cols$pwid]
                  ,RIGbiopax:::internal_seq_along_find_reps(gene_df_orig[,cols$pwid])
                  ,sep="_"
            )
        #make up controlled reaction ids in the same way
        #as control component ids
        gene_df_orig$rxn_id<-
            paste0(gene_df_orig$ctrl_id
                   ,"rxn")
        
        #split the pmid_author column lengthwise
        gene_df_orig<-
            gene_df_orig %>%
            split_cols_lengthen_df(colsToSplit = cols$pmid_author
                                   ,patternToSplit = "\\|"
            ) %>%
            #split pmid from author into sep columns
            split_col_widen_df(colToSplit = cols$pmid_author
                               ,split = ":"
                               ,newcolnames=c("pmid"
                                              ,"author"))

        #add pmids by merging original ctrl ids 
        #with "pmid" and additional numeric identifier to distinguish between different 
        #pmid/author combos that belong to the same entity
        gene_df_orig$pmid_id<-
            paste(gene_df_orig$ctrl_id
                  ,"pmid"
                  ,RIGbiopax:::internal_seq_along_find_reps(gene_df_orig$ctrl_id)
                  ,sep=""
            )
        
        #stack first gene df over second gene df
        gene_df<-
            data.frame(pwid=
                           gene_df_orig[,cols$pwid]
                       ,srcpwname=
                           gene_df_orig[,cols$srcpwname]
                       ,ctrl_id=
                           gene_df_orig$ctrl_id
                       ,ctrl_class=
                           gene_df_orig[,cols$ctrl_class]
                       ,ctrl_type=
                           gene_df_orig[,cols$ctrl_type]
                       ,rxn_id=
                           gene_df_orig$rxn_id
                       ,rxn_class=
                           gene_df_orig[,cols$rxn_class]
                       
                       ,ctrl_prop=
                           gene_df_orig$first_ctrl_prop
                       ,dispname=
                           gene_df_orig[,cols$first_dispname]
                       ,name=
                           gene_df_orig[,cols$first_name]
                       ,genesym=
                           gene_df_orig[,cols$first_genesym]
                       ,geneid=
                           gene_df_orig[,cols$first_geneid]
                       ,left_right=
                           gene_df_orig$first_left_right
                       ,pmid_id=
                           gene_df_orig$pmid_id
                       ,pmid=
                           gene_df_orig$pmid
                       ,author=
                           gene_df_orig$author
            ) %>%
            rbind.data.frame(data.frame(pwid=
                                            gene_df_orig[,cols$pwid]
                                        ,srcpwname=
                                            gene_df_orig[,cols$srcpwname]
                                        ,ctrl_id=
                                            gene_df_orig$ctrl_id
                                        ,ctrl_class=
                                            gene_df_orig[,cols$ctrl_class]
                                        ,ctrl_type=
                                            gene_df_orig[,cols$ctrl_type]
                                        ,rxn_id=
                                            gene_df_orig$rxn_id
                                        ,rxn_class=
                                            gene_df_orig[,cols$rxn_class]
                                        
                                        ,ctrl_prop=
                                            gene_df_orig$second_ctrl_prop
                                        ,dispname=
                                            gene_df_orig[,cols$second_dispname]
                                        ,name=
                                            gene_df_orig[,cols$second_name]
                                        ,genesym=
                                            gene_df_orig[,cols$second_genesym]
                                        ,geneid=
                                            gene_df_orig[,cols$second_geneid]
                                        ,left_right=
                                            gene_df_orig[,cols$second_left_right]
                                        ,pmid_id=
                                            gene_df_orig$pmid_id
                                        ,pmid=
                                            gene_df_orig$pmid
                                        ,author=
                                            gene_df_orig$author
                                        ))
        
        #split the name/geneid/etc cols by pipe
        gene_df<-
            gene_df %>%
            split_cols_lengthen_df(colsToSplit=c("dispname"
                                                 ,"name"
                                                 ,"genesym"
                                                 ,"geneid")
                                   ,patternToSplit = "\\|"
                                   ,at_once=TRUE
                                   )
        
        #which are not control
        noctrl<-
            is.na(gene_df$ctrl_class)
        #which are neither control nor reaction
        noctrlrxn<-
            is.na(gene_df$ctrl_class) &
            is.na(gene_df$rxn_class)
        
        
        #make individual entity (protein) ids
        #by merging name, displayname, and id
        #this way they'll be unique
        gene_df$entity_id<-
            gene_df$geneid
        #make up component ids
        #first, populate with control ids and classes
        gene_df$comp_id<-
            gene_df$ctrl_id
        gene_df$comp_class<-
            gene_df$ctrl_class
        #then, those which don't have control class, populate with rxn ids
        gene_df$comp_id[noctrl]<-
            gene_df$rxn_id[noctrl]
        gene_df$comp_class[noctrl]<-
            gene_df$rxn_class[noctrl]
        #finally, those which don't have ctrl class or rxn class, populate with entity ids
        #classes are protein
        gene_df$comp_id[noctrlrxn]<-
            gene_df$entity_id[noctrlrxn]
        gene_df$comp_class[noctrlrxn]<-
            "Protein"
        
        #add entityref and xref ids
        gene_df$entref_id<-
            paste0(gene_df$entity_id
                   ,"x")
        gene_df$xref_id<-
            paste0(gene_df$entity_id
                   ,"xx")
        
        ################################################################################
        #declare list to store all biopax-style data table
        dTable_list<-
            list()
        
        #pathway names
        dTable_list$dt_pw_names<-
            data.frame(class="Pathway"
                       ,id=gene_df$pwid
                       ,property="name"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df$srcpwname
            )
        
        #fill pathway components
        dTable_list$dt_comp<-
            data.frame(class="Pathway"
                       ,id=gene_df$pwid
                       ,property="pathwayComponent"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=gene_df$comp_id
                       ,property_value=""
            )
        #add components' pmids and author info
        dTable_list$dt_pmid<-
            data.frame(class=gene_df$comp_class
                       ,id=gene_df$comp_id
                       ,property="xref"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=gene_df$pmid_id
                       ,property_value=""
            ) 
        dTable_list$dt_pmid_db<-
            data.frame(class="PublicationXref"
                       ,id=gene_df$pmid_id
                       ,property="db"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value="Pubmed"
            ) 
        dTable_list$dt_pmid_id<-
            data.frame(class="PublicationXref"
                       ,id=gene_df$pmid_id
                       ,property="id"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df$pmid
            ) 
        dTable_list$dt_pmid_auth<-
            data.frame(class="PublicationXref"
                       ,id=gene_df$pmid_id
                       ,property="author"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df$author
            ) 
        
        #add control components' references
        #controllers
        dTable_list$dt_ctrl_r<-
            data.frame(class=gene_df$ctrl_class[gene_df$ctrl_prop=="controller"]
                       ,id=gene_df$ctrl_id[gene_df$ctrl_prop=="controller"]
                       ,property="controller"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=gene_df$entity_id[gene_df$ctrl_prop=="controller"]
                       ,property_value=""
            ) 
        
        #controlleds
        dTable_list$dt_ctrl_d<-
            data.frame(class=gene_df$ctrl_class[gene_df$ctrl_prop=="controlled"]
                       ,id=gene_df$ctrl_id[gene_df$ctrl_prop=="controlled"]
                       ,property="controlled"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=gene_df$rxn_id[gene_df$ctrl_prop=="controlled"]
                       ,property_value=""
            ) 
        
        #control types
        dTable_list$dt_ctrl_t<-
            data.frame(class=gene_df$ctrl_class
                       ,id=gene_df$ctrl_id
                       ,property="controlType"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df$ctrl_type
            ) 
        
        #add reactions refs
        dTable_list$dt_rxn<-
            data.frame(class=gene_df$rxn_class
                       ,id=gene_df$rxn_id
                       ,property=gene_df$left_right
                       ,property_attr="rdf:resource"
                       ,property_attr_value=gene_df$entity_id
                       ,property_value=""
            ) 
        
        #protein display names
        dTable_list$dt_prot_dispname<-
            data.frame(class="Protein"
                       ,id=gene_df$entity_id
                       ,property="displayName"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df$dispname
            )
        #protein names
        dTable_list$dt_prot_name<-
            data.frame(class="Protein"
                       ,id=gene_df$entity_id
                       ,property="name"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df$name
            )
        
        #protein entity references
        dTable_list$dt_prot_refs<-
            data.frame(class="Protein"
                       ,id=gene_df$entity_id
                       ,property="entityReference"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=gene_df$entref_id
                       ,property_value=""
            )
        
        #protein gene symbol
        dTable_list$dt_prot_refs_sym<-
            data.frame(class="ProteinReference"
                       ,id=gene_df$entref_id
                       ,property="name"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df$genesym
            )
        
        #protein entity references to the gene 
        dTable_list$dt_prot_refs_x<-
            data.frame(class="ProteinReference"
                       ,id=gene_df$entref_id
                       ,property="xref"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=gene_df$xref_id
                       ,property_value=""
            )
        
        #protein id db
        dTable_list$dt_prot_refs_x_db<-
            data.frame(class="RelationshipXref"
                       ,id=gene_df$xref_id
                       ,property="db"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value="entrezgene"
            )
        
        #protein gene id
        dTable_list$dt_prot_refs_x_id<-
            data.frame(class="RelationshipXref"
                       ,id=gene_df$xref_id
                       ,property="id"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value=gene_df$geneid
            )
        
        dTable<-
            do.call(rbind
                    ,dTable_list) %>%
            unique %>%
            .[complete.cases(.),]
        
        biopax<-
            biopax_from_dt(dTable
                           ,filename=filename)
        return(biopax)
        
    }



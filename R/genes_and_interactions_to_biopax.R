genes_and_interactions_to_biopax<-
  function(gene_df=NULL
           ,filename=NULL
           ,cols=list(pwid=
                        "toxdb.Pathway.ID"
                      ,srcpwname=
                        "source.Pathway.Name"
                      
                      ,first_class=
                        NULL
                      ,first_dispname=
                        "first_prot_sym"
                      ,first_name=
                        "first_name"
                      ,first_genesym=
                        "first_gene_sym"
                      ,first_geneid=
                        "first_gene_id"
                      ,first_species=
                        NULL
                      ,first_taxid=
                        NULL
                      
                      ,second_class=
                        NULL
                      ,second_dispname=
                        "second_prot_sym"	
                      ,second_name=
                        "second_name"
                      ,second_genesym=
                        "second_gene_sym"
                      ,second_geneid=
                        "second_gene_id"
                      ,second_species=
                        NULL
                      ,second_taxid=
                        NULL
                      
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
    #' @export
    #' @title
    #' Turn DataFrame of Interactions into BioPAX Object
    #' @description 
    #' Convert a dataframe containing a set of interactions into a BioPAX object.
    #' @param gene_df Dataframe containing a set of genes.
    #' @param filename (optional) Output filename.
    #' If provided, a resultant BioPAX object will also be written to file.
    #' @param cols Names of the pertaining expected columns in \code{gene_df}.
    
    #' @author 
    #' Ivan Grishagin
    
    #mandatory columns in the final df
    mandatory_cols<-
      c("pwid"
        ,"srcpwname"
        ,"ctrl_id"
        ,"ctrl_class"
        ,"ctrl_type"
        ,"rxn_id"
        ,"rxn_class"
        ,"ctrl_prop"
        ,"entity_class"
        ,"dispname"
        ,"name"
        ,"genesym"
        ,"geneid"
        ,"left_right"
        ,"pmid"
        ,"author"
        ,"species"
        ,"taxid")
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
                 
                 ,entity_class=
                   gene_df_orig[,cols$first_class]
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
                 ,pmid=
                   gene_df_orig$pmid
                 ,author=
                   gene_df_orig$author
                 
                 ,species=
                   gene_df_orig[,cols$first_species]
                 ,taxid=
                   gene_df_orig[,cols$first_taxid]
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
                                  
                                  ,entity_class=
                                    gene_df_orig[,cols$second_class]
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
                                  ,pmid=
                                    gene_df_orig$pmid
                                  ,author=
                                    gene_df_orig$author
                                  
                                  ,species=
                                    gene_df_orig[,cols$second_species]
                                  ,taxid=
                                    gene_df_orig[,cols$second_taxid]
      ))
    
    #fill columns of 0 length with NA values
    #0-length columns originate from assigining NULL-value columns
    for (cname in mandatory_cols) {
      if(length(gene_df[[cname]])<1){
        gene_df[[cname]]<-NA
      }
    }
    #make pmid ids by equating them to pmid/author combination
    #for those which don't have either -- fill NA
    gene_df$pmid_id<-
      paste0("pmid_"
             ,gene_df$pmid
      )
    gene_df$pmid_id[is.na(gene_df$pmid)]<-
      NA
    
    
    #replace all missing entity classes with "Protein"
    gene_df$entity_class[is.na(gene_df$entity_class)]<-
      "Protein"
    gene_df$entref_class<-
      paste0(gene_df$entity_class
             ,"Reference")
    
    #split the name/geneid/etc cols by pipe
    gene_df<-
      gene_df %>%
      split_cols_lengthen_df(colsToSplit=c("dispname"
                                           ,"name"
                                           ,"genesym"
                                           ,"geneid"
                                           ,"species"
                                           ,"taxid")
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
    #by displayname, and id
    #this way they'll be unique
    #just in case, if both disp name and gene id are NA
    #return NA
    gene_df$entity_id<-
      paste(gene_df$dispname
            ,gene_df$geneid
            ,sep="_")
    gene_df$entity_id[is.na(gene_df$dispname) &
                        is.na(gene_df$geneid)]<-
      NA
    
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
    #finally, those which don't have ctrl class or rxn class
    #populate with entity ids classes are protein
    gene_df$comp_id[noctrlrxn]<-
      gene_df$entity_id[noctrlrxn]
    gene_df$comp_class[noctrlrxn]<-
      gene_df$entity_class[noctrlrxn]
    
    #add entityref and xref ids
    #those, which do not have gene ids
    #make into NA values
    #so that there won't be links leading to nothing
    gene_df$entref_id<-
      paste0("protref_"
             ,gene_df$geneid
      )
    #if this is used, there will be multiple
    #xrefs for one entref -- the following 
    #extraction/remaking of biopax function does not cope with it
    # paste0("protref_"
    #        ,gene_df$genesym
    #        )
    
    gene_df$entref_id[is.na(gene_df$genesym)]<-
      NA
    gene_df$xref_id<-
      paste0("xref_"
             ,gene_df$geneid
      )
    gene_df$xref_id[is.na(gene_df$geneid)]<-
      NA
    
    #add organism and taxonomy ids
    gene_df$org_id<-
      paste0("biosrc_"
             ,gene_df$taxid)  
    gene_df$org_id[is.na(gene_df$species)]<-
      NA
    gene_df$taxid_id<-
      paste0("taxon_"
             ,gene_df$taxid) 
    gene_df$taxid_id[is.na(gene_df$taxid)]<-
      NA
    
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
      data.frame(class=gene_df$entity_class
                 ,id=gene_df$entity_id
                 ,property="displayName"
                 ,property_attr="rdf:datatype"
                 ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                 ,property_value=gene_df$dispname
      )
    #protein names
    dTable_list$dt_prot_name<-
      data.frame(class=gene_df$entity_class
                 ,id=gene_df$entity_id
                 ,property="name"
                 ,property_attr="rdf:datatype"
                 ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                 ,property_value=gene_df$name
      )
    
    #protein entity references
    dTable_list$dt_prot_refs<-
      data.frame(class=gene_df$entity_class
                 ,id=gene_df$entity_id
                 ,property="entityReference"
                 ,property_attr="rdf:resource"
                 ,property_attr_value=gene_df$entref_id
                 ,property_value=""
      )
    
    #protein gene symbol
    dTable_list$dt_prot_refs_sym<-
      data.frame(class=gene_df$entref_class
                 ,id=gene_df$entref_id
                 ,property="name"
                 ,property_attr="rdf:datatype"
                 ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                 ,property_value=gene_df$genesym
      )
    
    #protein entity references to the gene 
    dTable_list$dt_prot_refs_x<-
      data.frame(class=gene_df$entref_class
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
    
    #protein entity references to the gene 
    dTable_list$dt_org<-
      data.frame(class="ProteinReference"
                 ,id=gene_df$entref_id
                 ,property="organism"
                 ,property_attr="rdf:resource"
                 ,property_attr_value=gene_df$org_id
                 ,property_value=""
      )
    
    dTable_list$dt_org_name<-
      data.frame(class="BioSource"
                 ,id=gene_df$org_id
                 ,property="dispName"
                 ,property_attr="rdf:datatype"
                 ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                 ,property_value=gene_df$species
      )
    dTable_list$dt_org_tax_ref<-
      data.frame(class="BioSource"
                 ,id=gene_df$org_id
                 ,property="taxonXref"
                 ,property_attr="rdf:resource"
                 ,property_attr_value=gene_df$taxid_id
                 ,property_value=""
      )
    
    #protein id db
    dTable_list$dt_org_tax_db<-
      data.frame(class="UnificationXref"
                 ,id=gene_df$taxid_id
                 ,property="db"
                 ,property_attr="rdf:datatype"
                 ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                 ,property_value="entrezgene"
      )
    
    #protein gene id
    dTable_list$dt_org_tax_id<-
      data.frame(class="UnificationXref"
                 ,id=gene_df$taxid_id
                 ,property="id"
                 ,property_attr="rdf:datatype"
                 ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                 ,property_value=gene_df$taxid
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



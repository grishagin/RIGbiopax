add_symbols_entrezids_mult_keytypes<-
  function(dFrame
           ,col_names=
             c(id_col = "biopax.Gene.ID"
               ,type_col = "biopax.Gene.ID.Type"
               ,entrez_col = "ENTREZID"
               ,symbol_col = "biopax.Gene.Symbol"
               ,pathway_col = "biopax.Pathway.ID"
               ,comp_col = "biopax.Component.ID")
           ,keytypes=c("entrezgene"
                       ,"mim"
                       ,"uniprot"
                       ,"unigene"
                       ,"ensemblgene"
                       ,"hprd"
                       ,"hgnc")
           ,filter_keytypes=TRUE){
    #' @export
    #' @title 
    #' Annotate DataFrame with Gene ID and Gene Symbol
    #' @description 
    #' Adds a Gene ID and Symbol columns to a dataframe with a protein/gene 
    #' ids in multiple other formats.
    #' @param dFrame Dataframe or data table to process.
    #' @param col_names A named vector of actual column names 
    #' in the supplied \code{dFrame}. 
    #' @param keytypes A format of IDs to be converted. 
    #' @param filter_keytypes Logical. Return only rows with desired 
    #' format types? Defaults to \code{TRUE}.
    #' 
    #' @author 
    #' Ivan Grishagin
    
    #ensure it's a data table
    dFrame<-
      dFrame %>% 
      as.data.table
    
    #replace the key column names 
    colnames(dFrame)<-
      colnames(dFrame) %>% 
      mapvalues(from=col_names
                ,to=names(col_names)
                ,warn_missing = FALSE)
    
    #if entrez or symbol cols are not present -- add them
    if(is.null(dFrame$entrez_col)){
      dFrame$entrez_col<-NA_integer_
    }
    if(is.null(dFrame$symbol_col)){
      dFrame$symbol_col<-NA
    }
    
    #ensure correct types of Gene Symbol and Entrezid columns
    dFrame$symbol_col<-
      as.character(dFrame$symbol_col)
    
    dFrame$entrez_col<-
      as.integer(dFrame$entrez_col)
    
    
    
    #if there are multiple comma|semicolon|pipe-separated values,
    #split them and expand the df
    dFrame<-
      dFrame %>%
      split_cols_lengthen_df(colsToSplit="id_col"
                             ,patternToSplit = ",|;|\\|"
                             ,at_once=TRUE) %>% 
      as.data.table
    
    #annotate data table iteratively
    for(KEYTYPE in keytypes){
      message("Querying ", KEYTYPE," keytype using mygene...")
      
      dFrame<-
        dFrame %>%
        add_symbols_entrezids_single_keytype(col_names = 
                                               c(id_col = "id_col"
                                                 ,type_col = "type_col"
                                                 ,entrez_col = "entrez_col"
                                                 ,symbol_col = "symbol_col"
                                                 ,pathway_col = "pathway_col"
                                                 ,comp_col = "comp_col")
                                             ,KEYTYPE = KEYTYPE)
    }
    
    #if multiple pipe-separated values were added when annotating,
    #split them and expand the df
    dFrame<-
      dFrame %>%
      split_cols_lengthen_df(colsToSplit=c("symbol_col"
                                           ,"entrez_col")
                             ,patternToSplit = "\\|"
                             ,at_once=TRUE) %>% 
      as.data.table
    
    
    if(filter_keytypes){
      dFrame<-
        dFrame[type_col %chin% keytypes]
    }
    
    #change colnames back
    colnames(dFrame)<-
      colnames(dFrame) %>% 
      mapvalues(from=names(col_names)
                ,to=col_names
                ,warn_missing = FALSE)
    
    #return the annotated df
    return(dFrame)
    
  }
############################################################################################
# internal_add_multiple_symbols_entrezids <-
#     function(df_pw_proteins=NULL
#              ,keytypes=c("entrezgene"
#                          ,"mim"
#                          ,"uniprot"
#                          ,"unigene"
#                          ,"ensemblgene"
#                          ,"hprd"
#                          ,"hgnc")
#              ,filter_keytypes=TRUE){
#         
#         require(plyr)
#         require(dplyr)
#         
#         if (is.null(df_pw_proteins)){
#             stop("add_symbols_entrezids: df_pw_proteins is missing!")
#         }
#         
#         #replace db ids name with proper annotation ids
#         df_pw_proteins$biopax.Gene.ID.Type<-
#             tolower(df_pw_proteins$biopax.Gene.ID.Type)
#         df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "ll"]<-"entrezgene"
#         df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "entrezgene"]<-"entrezgene"
#         df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "entrez gene"]<-"entrezgene"
#         df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "uniprot"]<-"uniprot"
#         df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "unigene"]<-"unigene"
#         df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "ensembl"]<-"ensemblgene"
#         df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "omim"]<-"mim"
#         df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "hprd"]<-"hprd"
#         df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "hgnc"]<-"hgnc"
#         
#         
#         #ensure correct types of Gene Symbol and Entrezid columns
#         df_pw_proteins$biopax.Gene.Symbol<-
#             as.character(df_pw_proteins$biopax.Gene.Symbol)
#         
#         df_pw_proteins$ENTREZID<-
#             as.integer(df_pw_proteins$ENTREZID)
#         
#         #if there are multiple comma|semicolon|pipe-separated values,
#         #split them and expand the df
#         df_pw_proteins<-
#             df_pw_proteins %>%
#             split_cols_lengthen_df(colsToSplit="biopax.Gene.ID"
#                                    ,patternToSplit = ",|;|\\|"
#                                    ,at_once=TRUE)
#         
#         #annotate data table iteratively
#         for(KEYTYPE in keytypes){
#             df_pw_proteins<-
#                 df_pw_proteins %>%
#                 internal_add_symbols_entrezids_mygene(KEYTYPE)
#         }
#         
#         #if multiple pipe-separated values were added when annotating,
#         #split them and expand the df
#         df_pw_proteins<-
#             df_pw_proteins %>%
#             split_cols_lengthen_df(colsToSplit=c("biopax.Gene.Symbol"
#                                                  ,"ENTREZID")
#                                    ,patternToSplit = "\\|"
#                                    ,at_once=TRUE)
#         
#         
#         if(filter_keytypes){
#             df_pw_proteins<-
#                 df_pw_proteins %>%
#                 filter(biopax.Gene.ID.Type %chin% keytypes)
#         }
#         
#         #return the annotated df
#         return(df_pw_proteins)
#         
#     }

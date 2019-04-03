add_symbols_entrezids_single_keytype<-
  function(dFrame
           ,col_names=
             c(id_col = "id_col"
               ,type_col = "type_col"
               ,entrez_col = "entrez_col"
               ,symbol_col = "symbol_col"
               ,pathway_col = ""
               ,comp_col = "")
           ,KEYTYPE=NULL
           ,species="human"){
    #' @export
    #' @title 
    #' Annotate DataFrame with Gene ID and Gene Symbol
    #' @description 
    #' Adds a Gene ID and Symbol columns to a dataframe with a protein/gene ids 
    #' in one other format (see \code{KEYTYPE} argument).
    #' @param dFrame Dataframe or data table to process.
    #' @param col_names A named vector of actual column names 
    #' in the supplied \code{dFrame}. 
    #' @param KEYTYPE A format of IDs to be converted. 
    #' @param species Species to look for. Default is "human".
    #' 
    #' @author 
    #' Ivan Grishagin
    
    #ensure it's a data table
    dFrame<-
      dFrame %>% 
      as.data.table
    
    #set default values for query results
    queryResults<-
      NULL
    #replace the key column names 
    colnames(dFrame)<-
      colnames(dFrame) %>% 
      mapvalues(from=col_names
                ,to=names(col_names)
                ,warn_missing = FALSE)
    
    #tolower all type column values
    dFrame$type_col<-
      dFrame$type_col %>% 
      tolower
    
    #fix certain keytypes since biopax has a variety of db names 
    to_replace<-
      switch(KEYTYPE
             ,entrezgene = c("ll"
                             ,"entrez gene")
             ,mim = "omim"
             ,ensemblgene = "ensembl")
    
    dFrame[type_col %in% to_replace]$type_col<-
      KEYTYPE
    
    #select correct rows
    KEYTYPE_rows<-
      which(dFrame$type_col %chin% KEYTYPE)
    
    #if such rows are available, convert said type to ENTREZID or Symbol
    if(length(KEYTYPE_rows)>0){
      
      #lookup pertaining gene ids to convert
      queryInput<-
        dFrame$id_col[KEYTYPE_rows]
      
      if(KEYTYPE == "entrezgene"){
        #if querying entrezgenes -- use a different function
        #this way, it will return results regardless of species
        #otherwise, substantial modification of the script is required
        #to accommodate for different species
        
        queryResults<-
          try(getGenes(unique(queryInput)
                       ,fields=c("entrezgene"
                                 ,"symbol")
                       ,return.as="DataFrame"))
      } else {
        #query the input values
        queryResults<-
          try(queryMany(unique(queryInput)
                        ,scopes=KEYTYPE
                        ,fields=c("entrezgene"
                                  ,"symbol")
                        ,return.as="DataFrame"
                        ,species=species))
      }
      
      
      if(!("entrezgene" %chin% colnames(queryResults))){
        message("Unsuccessful query for ", KEYTYPE,".")
        return(dFrame)
      } else {
        message("Successful query for ", KEYTYPE,"!")
      }
      
      queryResults<-
        queryResults %>%
        merge_cols_shorten_df(colKey = "query"
                              ,colsToMerge = c("entrezgene"
                                               ,"symbol")
                              ,patternToMerge="|")
      
      #fill the entrez gene and symbol values
      dFrame[KEYTYPE_rows]$entrez_col<-
        queryResults$entrezgene[match(queryInput
                                      ,queryResults$query)] %>% 
        as.integer
      dFrame[KEYTYPE_rows]$symbol_col<-
        queryResults$symbol[match(queryInput
                                  ,queryResults$query)]
      
      
    }
    
    ############################## special netpath case ############################## 
    ################################################################################## 
    #for netpath, try to extract symbols from element ids
    if(all(c("pathway_col","comp_col") %in% colnames(dFrame))){
      if(length(grep("netpath"
                     ,unique(dFrame$pathway_col)
                     ,ignore.case = TRUE))>0){
        
        #find rows with NA values for a particular keytype
        KEYTYPE_NA_rows<-
          which(dFrame$type_col %chin% KEYTYPE &
                  is.na(dFrame$id_col))
        #set default values for query results
        queryResults<-
          NULL
        
        if(length(KEYTYPE_NA_rows)>0){
          #for those rows, extract gene symbols from the component ID column
          symbols<-
            dFrame$comp_col[KEYTYPE_NA_rows] %>%
            strsplit(split="_") %>%
            sapply("[[",2)
          
          #try to get entrez id for each of the symbols
          while(is.null(queryResults)){
            #sometimes curl error is thrown, so try until it budges
            queryResults<-
              try(queryMany(unique(symbols)
                            ,scopes="symbol"
                            ,fields="entrezgene"
                            ,return.as="DataFrame"
                            ,species=species))
            # %>% 
            #     data.frame 
            
            if(class(queryResults)=="try-error"){
              queryResults<-NULL
            }
          }
          
          
          if(!("entrezgene" %chin% colnames(queryResults))){
            message("Unsuccessful query for symbols: "
                    ,unique(symbols))
            return(dFrame)
          }
          
          queryResults<-
            queryResults %>%
            merge_cols_shorten_df(colKey = "query"
                                  ,colsToMerge = "entrezgene"
                                  ,patternToMerge="|")
          
          dFrame[KEYTYPE_rows]$entrez_col<-
            queryResults$entrezgene[match(symbols
                                          ,queryResults$query)]
          dFrame[KEYTYPE_rows]$symbol_col<-
            queryResults$query[match(symbols
                                     ,queryResults$query)]
        }
      }
      
    }
    
    #change colnames back
    colnames(dFrame)<-
      colnames(dFrame) %>% 
      mapvalues(from=names(col_names)
                ,to=col_names
                ,warn_missing = FALSE)
    
    return(dFrame)
  }
###############################################################################################
# internal_add_symbols_entrezids_mygene<-
#     function(df_pw_proteins=NULL
#              ,KEYTYPE=NULL
#              ,species="human"){
#         require(plyr)
#         require(dplyr)
#         require(mygene)
#         
#         if (is.null(df_pw_proteins) |
#             is.null(KEYTYPE)){
#             stop("add_symbols_entrezids: df_pw_proteins or KEYTYPE is missing!")
#         }
#         
#         #select correct rows
#         KEYTYPE_rows<-
#             which(df_pw_proteins$biopax.Gene.ID.Type %chin% KEYTYPE)
#         
#         #if such rows are available, convert said type to ENTREZID or Symbol
#         if(length(KEYTYPE_rows)>0){
#             
#             #lookup pertaining gene ids to convert
#             queryInput<-
#                 df_pw_proteins$biopax.Gene.ID[KEYTYPE_rows]
#             
#             if(KEYTYPE=="entrezgene"){
#                 #if querying entrezgenes -- use a different function
#                 #this way, it will return results regardless of species
#                 #otherwise, substantial modification of the script is required
#                 #to accommodate for different species
#                 queryResults<-
#                     getGenes(unique(queryInput)
#                              ,fields=c("entrezgene","symbol")
#                              ,return.as="DataFrame") %>% 
#                     data.frame
#             } else {
#                 #query the input values
#                 queryResults<-
#                     queryMany(unique(queryInput)
#                               ,scopes=KEYTYPE
#                               ,fields=c("entrezgene","symbol")
#                               ,return.as="DataFrame"
#                               ,species=species) %>% 
#                     data.frame
#             }
#             
#             
#             if(!("entrezgene" %chin% colnames(queryResults))){
#                 message("Unsuccessful query for ", KEYTYPE,".")
#                 return(df_pw_proteins)
#             }
#             
#             queryResults<-
#                 queryResults %>%
#                 merge_cols_shorten_df(colKey = "query"
#                                       ,colsToMerge = c("entrezgene","symbol")
#                                       ,patternToMerge="|")
#             
#             #fill the entrez gene and symbol values
#             df_pw_proteins[KEYTYPE_rows
#                            ,c("ENTREZID","biopax.Gene.Symbol")]<-
#                 queryResults[match(queryInput,queryResults$query),
#                              c("entrezgene","symbol")]
#             
#         }
#         
#         #for netpath
#         if(length(grep("netpath"
#                        ,unique(df_pw_proteins$biopax.Pathway.ID)
#                        ,ignore.case = TRUE))>0){
#             
#             #find rows with NA values for a particular keytype
#             KEYTYPE_NA_rows<-
#                 which(df_pw_proteins$biopax.Gene.ID.Type %chin% KEYTYPE &
#                           is.na(df_pw_proteins$biopax.Gene.ID))
#             
#             if(length(KEYTYPE_NA_rows)>0){
#                 #for those rows, extract gene symbols from the component ID column
#                 symbols<-
#                     df_pw_proteins$biopax.Component.ID[KEYTYPE_NA_rows] %>%
#                     strsplit(split="_") %>%
#                     sapply("[[",2)
#                 
#                 #try to get entrez id for each of the symbols
#                 queryResults<-
#                     queryMany(unique(symbols)
#                               ,scopes="symbol"
#                               ,fields="entrezgene"
#                               ,return.as="DataFrame"
#                               ,species=species) %>% 
#                     data.frame 
#                 
#                 if(!("entrezgene" %chin% colnames(queryResults))){
#                     message("Unsuccessful query for symbols: ",unique(symbols))
#                     return(df_pw_proteins)
#                 }
#                 
#                 queryResults<-
#                     queryResults %>%
#                     merge_cols_shorten_df(colKey = "query"
#                                           ,colsToMerge = "entrezgene"
#                                           ,patternToMerge="|")
#                 
#                 df_pw_proteins[KEYTYPE_NA_rows
#                                ,c("ENTREZID","biopax.Gene.Symbol")]<-
#                     queryResults[match(symbols,queryResults$query),
#                                  c("entrezgene","query")]
#             }
#         }
#         
#         
#         return(df_pw_proteins)
#     }

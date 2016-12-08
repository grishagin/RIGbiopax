add.symbols.entrezids.mygene <-
function(df_pw_proteins=NULL
         ,KEYTYPE=NULL
         ,species="human"
    ){
        require(plyr)
        require(dplyr)
        require(mygene)
        
        if (is.null(df_pw_proteins) |
            is.null(KEYTYPE)){
            stop("add_symbols_entrezids: df_pw_proteins or KEYTYPE is missing!")
        }
        
        #select correct rows
        KEYTYPE_rows<-
            which(df_pw_proteins$biopax.Gene.ID.Type %chin% KEYTYPE)
        
        #if such rows are available, convert said type to ENTREZID or Symbol
        if(length(KEYTYPE_rows)>0){
            
            #lookup pertaining gene ids to convert
            queryInput<-
                df_pw_proteins$biopax.Gene.ID[KEYTYPE_rows]
            
            if(KEYTYPE=="entrezgene"){
                #if querying entrezgenes -- use a different function
                #this way, it will return results regardless of species
                #otherwise, substantial modification of the scrip is required
                #to accommodate for different species
                queryResults<-
                    getGenes(unique(queryInput)
                             ,fields=c("entrezgene","symbol")
                             ,return.as="DataFrame") %>% 
                    data.frame
            } else {
                #query the input values
                queryResults<-
                    queryMany(unique(queryInput)
                              ,scopes=KEYTYPE
                              ,fields=c("entrezgene","symbol")
                              ,return.as="DataFrame"
                              ,species=species) %>% 
                    data.frame
            }

            
            if(!("entrezgene" %chin% colnames(queryResults))){
                message("Unsuccessful query for ", KEYTYPE,".")
                return(df_pw_proteins)
            }
            
            queryResults<-
                queryResults %>%
                shrink.df.via.merge.col(colKey = "query"
                                        ,colToMerge = c("entrezgene","symbol")
                                        ,patternToMerge="|")
            
            #fill the entrez gene and symbol values
            df_pw_proteins[KEYTYPE_rows
                           ,c("ENTREZID","biopax.Gene.Symbol")]<-
                queryResults[match(queryInput,queryResults$query),
                             c("entrezgene","symbol")]
            
        }
        
        #for netpath
        if(length(grep("netpath"
                       ,unique(df_pw_proteins$biopax.Pathway.ID)
                       ,ignore.case = TRUE))>0){
            
            #find rows with NA values for a particular keytype
            KEYTYPE_NA_rows<-
                which(df_pw_proteins$biopax.Gene.ID.Type %chin% KEYTYPE &
                          is.na(df_pw_proteins$biopax.Gene.ID))
            
            if(length(KEYTYPE_NA_rows)>0){
                #for those rows, extract gene symbols from the component ID column
                symbols<-
                    df_pw_proteins$biopax.Component.ID[KEYTYPE_NA_rows] %>%
                    strsplit(split="_") %>%
                    sapply("[[",2)
                
                #try to get entrez id for each of the symbols
                queryResults<-
                    queryMany(unique(symbols)
                              ,scopes="symbol"
                              ,fields="entrezgene"
                              ,return.as="DataFrame"
                              ,species=species) %>% 
                    data.frame 
                
                if(!("entrezgene" %chin% colnames(queryResults))){
                    message("Unsuccessful query for symbols: ",unique(symbols))
                    return(df_pw_proteins)
                }
                
                queryResults<-
                    queryResults %>%
                    shrink.df.via.merge.col(colKey = "query"
                                            ,colToMerge = "entrezgene"
                                            ,patternToMerge="|")
                
                df_pw_proteins[KEYTYPE_NA_rows
                               ,c("ENTREZID","biopax.Gene.Symbol")]<-
                    queryResults[match(symbols,queryResults$query),
                                 c("entrezgene","query")]
            }
        }
        
        
        return(df_pw_proteins)
    }

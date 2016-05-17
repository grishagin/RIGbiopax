compare_pw_components <-
function(biopax_prot=NULL
             ,toxdb_prot=NULL
             ,dfrow=NULL
             ,output=c("all","match","biopax","toxdb")
             
    ){
        #prepare a dummy dataframe
        #to be used in case the supplied pathways don't have any components
        dummy_df<-
            adjust_columns(df_to_adjust = dfrow
                           ,df_template = biopax_prot)
        #get toxdb pw id and biopax pw name
        pw_toxdb_id=dfrow$toxdb.Pathway.ID
        pw_biopax_id=dfrow$biopax.Pathway.ID
        
        #check if either of pathway names was NA
        #if so, return the other dataframe
        if (is.na(pw_biopax_id) &
            is.na(pw_toxdb_id)){
            message("compare_pw_components: both pw_biopax_id and pw_toxdb_id are NA!")
            return(dummy_df)
        }
        if (is.na(pw_biopax_id)){
            toxdb_prot_sub<-
                toxdb_prot %>%
                filter(toxdb.Pathway.ID==pw_toxdb_id)
            if(nrow(biopax_prot_sub)<1){
                message("compare_pw_components: biopax_prot_sub is NA toxdb_prot_sub has 0 rows!")
                return(dummy_df)
            }
            return(toxdb_prot_sub)
        } else if (is.na(pw_toxdb_id)){
            biopax_prot_sub<-
                biopax_prot %>%
                filter(biopax.Pathway.ID==pw_biopax_id)
            if(nrow(biopax_prot_sub)<1){
                message("compare_pw_components: pw_toxdb_id is NA and biopax_prot_sub has 0 rows!")
                return(dummy_df)
            }
            return(biopax_prot_sub)
        }
        #get subsets of the dataframes
        #containing components of corresponding pathway in question
        #for both data table (toxdb, etc.) and biopax
        
        biopax_prot_sub<-
            biopax_prot %>%
            filter(biopax.Pathway.ID==pw_biopax_id & toxdb.Pathway.ID==pw_toxdb_id)
        
        #df_prot_sub pw_bio_name pw_df_name df_prot bio_prot
        #fill biopax credentials in toxdb subframe
        toxdb_prot_sub<-
            toxdb_prot %>%
            filter(toxdb.Pathway.ID==pw_toxdb_id) %>%
            mutate(biopax.Pathway.Name = as.character(dfrow$biopax.Pathway.Name)
                   ,biopax.Pathway.ID = as.character(dfrow$biopax.Pathway.ID)
                   ,pathway.Match.Status = as.character(dfrow$pathway.Match.Status))
        
        if(nrow(biopax_prot_sub)<1){
            message("compare_pw_components: biopax_prot_sub has 0 rows!")
            return(toxdb_prot_sub)
        } else if(nrow(toxdb_prot_sub)<1){
            message("compare_pw_components: toxdb_prot_sub has 0 rows!")
            return(biopax_prot_sub)
        } else if (nrow(biopax_prot_sub)<1 &
                   nrow(toxdb_prot_sub)<1){
            message("compare_pw_components: biopax_prot_sub and toxdb_prot_sub have 0 rows!")
            return(dummy_df)
        } 
        
        #match which gene ids from biopax and data table are present in both
        pres_toxdb.pres_biopax<-
            biopax_prot_sub %>%
            .[(.$ENTREZID %in% toxdb_prot_sub$toxdb.Gene.ID),]
        pres_toxdb.pres_biopax$toxdb.Gene.ID<-
            toxdb_prot_sub$toxdb.Gene.ID[match(pres_toxdb.pres_biopax$ENTREZID
                                               ,toxdb_prot_sub$toxdb.Gene.ID)]
        pres_toxdb.pres_biopax$toxdb.Gene.Symbol<-
            toxdb_prot_sub$toxdb.Gene.Symbol[match(pres_toxdb.pres_biopax$ENTREZID
                                                   ,toxdb_prot_sub$toxdb.Gene.ID)]
        
        #match which gene ids from biopax are missing from data table
        pres_biopax.miss_toxdb<-
            biopax_prot_sub %>%
            .[!(.$ENTREZID %in% toxdb_prot_sub$toxdb.Gene.ID),]
        
        #...and which rows from data table are missing from biopax
        pres_toxdb.miss_biopax<-
            toxdb_prot_sub %>%
            .[!(.$toxdb.Gene.ID %in% biopax_prot_sub$ENTREZID),]
        
        if(output=="all"){
            to_return<-
                rbind(pres_toxdb.pres_biopax
                      ,pres_biopax.miss_toxdb
                      ,pres_toxdb.miss_biopax)
            return(to_return)
        } else if (output=="match"){
            return(pres_toxdb.pres_biopax)
            
        } else if (output=="biopax") {
            return(pres_biopax.miss_toxdb)
            
        } else if (output=="toxdb") {
            return(pres_toxdb.miss_biopax)
        }
        
    }

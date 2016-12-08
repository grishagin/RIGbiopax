internal_rm_duplicated_biopax_pw<-
    function(big_df){
        
        require(plyr)
        require(dplyr)
        
        #add gene match indicators
        .genematch<-
            rep(NA,nrow(big_df))
        #match
        .genematch[which(big_df$ENTREZID==big_df$toxdb.Gene.ID)]<-
            "1_match"
        #only in toxdb list
        .genematch[which(is.na(big_df$biopax.Gene.ID.Type) & !is.na(big_df$toxdb.Gene.ID))]<-
            "2_toxdb"
        #only in biopax
        .genematch[which(!is.na(big_df$biopax.Gene.ID.Type) & is.na(big_df$toxdb.Gene.ID))]<-
            "3_biopax"
        #neither in biopax nor in the toxdb list
        .genematch[which(is.na(big_df$biopax.Gene.ID.Type) & is.na(big_df$toxdb.Gene.ID))]<-
            "4_genes_not_found"
        #add column to the big df with genes
        big_df <-
            big_df%>%
            mutate(.genematch=.genematch) 

        #filter out individual toxdb/biopax id pairs
        #and match counts
        comp_df<-
            big_df %>%
            filter(!is.na(toxdb.Pathway.ID)) %>%
            group_by(toxdb.Pathway.ID
                     ,biopax.Pathway.ID
                     ,.genematch) %>%
            dplyr::summarise(count=n()) %>%
            spread(key=.genematch
                   ,value=count
                   ,fill = 0) %>%
            mutate(sum=`1_match`+`2_toxdb`+`3_biopax`)
        
        dupl_txid<-
            comp_df$toxdb.Pathway.ID[duplicated(comp_df$toxdb.Pathway.ID)]

        #for each duplicate txid, find biopax ids to be removed
        #because there has to be only one!
        bad_bpids<-
            dupl_txid %>%
            lapply(function(txid){
                #filter by tox id
                #and arrange by max match, then max sum
                #take all but the first bp ids
                bad_bpid<-
                    comp_df %>%
                    filter(toxdb.Pathway.ID %in% txid) %>%
                    # dplyr::arrange(desc(`1_match`)
                    #                ,desc(sum)
                    #                ) %>%
                    .$biopax.Pathway.ID %>%
                    .[-1]
                return(bad_bpid)
            }) %>%
            unlist
        
        #clean big df
        big_df<-
            big_df %>%
            filter(!biopax.Pathway.ID %in% bad_bpids)
        
        return(big_df)
    }
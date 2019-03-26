internal_analyse_toxdb_biopax_results_df<-
    function(dFrame
             ,label=NULL){

        #make a summary
        summary_df<-
            dFrame %>%
            #group by pathway id, then genematch
            dplyr::group_by(toxdb.Pathway.ID,.genematch) %>%
            dplyr::summarise(count=n()) %>%
            spread(key = .genematch
                   ,value = count) %>%
            dplyr::mutate(total_toxdb=sum(`1_match`
                                          ,`2_toxdb`
                                          ,na.rm=TRUE)
                          ,total_biopax=sum(`1_match`
                                            ,`3_biopax`
                                            ,na.rm=TRUE)
                          ,`matches_toxdb,%`=
                              round(100*`1_match`/total_toxdb,0) 
                          ,`misses_toxdb,%`=
                              round(100*`2_toxdb`/total_toxdb,0)
                          ,`matches_biopax,%`=
                              round(100*`1_match`/total_biopax,0)
                          ,`extra_biopax,%`=
                              round(100*`3_biopax`/total_biopax,0)
                          ) 
        #replace NA with 0
        summary_df[is.na(summary_df)]<-0
        
        #arrange by misses_toxdb
        summary_df<-
            summary_df %>%
            dplyr::arrange(`matches_toxdb,%`
                           ,`matches_biopax,%`
                           ,`extra_biopax,%`)
        
        #fix colnames
        colnames(summary_df)[2:4]<-
            c("matches"
              ,"misses_toxdb"
              ,"extra_biopax")
        
        #write to file
        write.xlsx(summary_df
                   ,paste0(Sys.Date()
                          ,"_toxdb_biopax_summary_analysis--"
                          ,label
                          ,"--genes.xlsx")
                   ,col.names = TRUE
                   ,row.names = FALSE)    
        
        return(summary_df)
           
    }







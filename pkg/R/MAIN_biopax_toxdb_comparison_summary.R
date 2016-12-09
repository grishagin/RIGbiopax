MAIN_biopax_toxdb_comparison_summary<-
    function(work_dir="D:/Dropbox/Rancho/NCATS/ToxDB/2016-12-08 RESULTS toxdb-vs-biopax genes comparison"
             ,output_dir="./summary"
             ,pw_matchup_file="../_source_files/pathways_matched_to_sources_v015.xlsx"
             ,pw_src_file="../_source_files/human_pathways_rbs_curated_20150615_with_source_info.xlsx"
             ,gene_df_file="./2016-12-08_ALL_SOURCES_toxdb_vs_biopax_genes.xlsx"
             ){
        
        #' @title
        #' Prepare BioPAX-ToxDB Comparison Summary
        #' @description 
        #' Prepares comprehensive summary of comparison of genes in BioPAX sources to ToxDB gene list.
        #' @param work_dir Work directory.
        #' @param output_dir Output directory (will be created if non-existent).
        #' @param pw_matchup_file File with pathway matchups between BioPAX sources and ToxDB list.
        #' @param pw_src_file ToxDB pathway spreadsheet with curation information.
        #' @param gene_df_file Spreadsheet with BioPAX-ToxDB gene comparison.
        
        #' @author 
        #' Ivan Grishagin
        
        require(RIGessentials)
        require(readxl)
        prepareSession(work_dir)
        
        dir.create(output_dir
                   ,showWarnings = FALSE)
        
        pw_match<-
            read_excel(path = pw_matchup_file
                       ,col_types = rep("text",11)
                       ,sheet = 1) 
        
        pw_src<-
            read_excel(path = pw_src_file
                       ,col_types = rep("text",4)
                       ,sheet = 1) 
        
        
        gene_df<-
            read_excel(path = gene_df_file
                       ,col_types = rep("text",14)
                       ,sheet = 1) %>%
            filter(!is.na(toxdb.Pathway.ID))
        
        #non-curated (non-altered) pathways
        pw_noncur<-
            pw_src %>%
            filter(Changed=="No") %>%
            .$Pathway
        
        #pathways with alternative source
        pw_altsrc<-
            pw_match %>%
            filter(!is.na(proposed.Source)) %>%
            .$toxdb.Pathway.ID
        
        #list of summary outputs
        summary_list<-list()
        
        #summary of original dataframe 
        summary_list$all<-
            gene_df %>%
            RIGbiopax:::internal_analyse_toxdb_biopax_results_df("all")
        
        #summary of dataframe with pathways only from original sources
        summary_list$origsrc<-
            gene_df %>%
            filter(!toxdb.Pathway.ID %in% pw_altsrc) %>%
            RIGbiopax:::internal_analyse_toxdb_biopax_results_df("origsrc")
        
        #summary of dataframe with pathways only from original sources+exact name matches
        summary_list$origsrc_exactmatch<-
            gene_df %>%
            filter(!toxdb.Pathway.ID %in% pw_altsrc) %>%
            filter(pathway.Match.Status=="Exact") %>%
            RIGbiopax:::internal_analyse_toxdb_biopax_results_df("origsrc_exactmatch")
        
        #summary of dataframe with pathways only from original sources+exact name matches+not changed
        summary_list$origsrc_exactmatch_noncur<-
            gene_df %>%
            filter(!toxdb.Pathway.ID %in% pw_altsrc) %>%
            filter(pathway.Match.Status=="Exact") %>%
            filter(toxdb.Pathway.ID %in% pw_noncur) %>%
            RIGbiopax:::internal_analyse_toxdb_biopax_results_df("origsrc_exactmatch_noncur")
        
        summary_cols<-
            c("all"
              ,"only orig src"
              ,"only orig src + exact matches"
              ,"only orig src + exact matches + not changed"
            )
        
        total_genes<-
            summary_list %>%
            lapply(FUN=function(dFrame){
                summ_vect<-
                    nrow(dFrame)
                names(summ_vect)<-"N_genes"
                summ_vect<-
                    c(summ_vect
                      ,dFrame[,c("matches"
                                 ,"misses_toxdb"
                                 ,"extra_biopax")] %>%
                          apply(MARGIN = 2
                                ,sum))
                return(summ_vect)
            }) %>%
            do.call(rbind,.) %>%
            cbind.data.frame(Subset=summary_cols
                             ,.) %>%
            dplyr::mutate(`misses_toxdb,%`=round(100*misses_toxdb/(misses_toxdb+
                                                                matches)
                                                 ,0))
        
        
        #quantiles for all subsets
        quantile_summary_df<-
            summary_list %>%
            lapply(RIGbiopax:::internal_prepare_quantile_summary) %>%
            do.call(rbind.data.frame,.) %>%
            cbind.data.frame(Subset=rep(summary_cols
                                        ,each=4)
                             ,.)
        
        #pathways with no genes in biopax
        pw_nogenes<-
            gene_df %>%
            dplyr::group_by(toxdb.Pathway.ID) %>%
            dplyr::summarise(na.len=sum(is.na(biopax.Component.ID))
                      ,all.len=length(biopax.Component.ID)
                      ,diff=all.len-na.len
                      ,Source=unique(Source)) %>%
            filter(diff==0) %>%
            dplyr::select(toxdb.Pathway.ID,Source) 
        
        
        write.xlsx(total_genes
                   ,paste0(Sys.Date()
                           ,"_toxdb_biopax_summary_analysis--total_genes.xlsx") %>% 
                       file.path(output_dir
                                 ,.)
                   ,col.names = TRUE
                   ,row.names = FALSE) 
        
        
        write.xlsx(quantile_summary_df
                   ,paste0(Sys.Date()
                           ,"_toxdb_biopax_summary_analysis--quantiles.xlsx") %>% 
                       file.path(output_dir
                                 ,.)
                   ,col.names = TRUE
                   ,row.names = FALSE) 
        
        
        write.xlsx(pw_nogenes
                   ,paste0(Sys.Date()
                           ,"_toxdb_biopax_summary_analysis--pw_no_genes.xlsx") %>% 
                       file.path(output_dir
                                 ,.)
                   ,col.names = TRUE
                   ,row.names = FALSE) 
        
        
        
    }
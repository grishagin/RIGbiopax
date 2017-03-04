MAIN_biopax_toxdb_comparison_summary<-
    function(output_dir="./summary"
             ,pw_matchup_file="default"
             ,pw_src_file="default"
             ,genes_df
             ){
        
        #' @title
        #' Prepare BioPAX-ToxDB Comparison Summary
        #' @description 
        #' Prepares comprehensive summary of comparison of genes in BioPAX sources to ToxDB gene list.
        #' @param output_dir Output directory (will be created if non-existent).
        #' @param pw_matchup_file File with pathway matchups between BioPAX sources and Inxight Pathways list.
        #' @param pw_src_file Inxight Pathways spreadsheet with pathway curation information.
        #' @param genes_df Either a dataframe with BioPAX-ToxDB gene comparison 
        #' or a filename of a spreadsheet with such dataframe.
        
        #' @author 
        #' Ivan Grishagin
        
        require(RIGessentials)
        require(readxl)
        
        dir.create(output_dir
                   ,showWarnings = FALSE)
        prepareSession(output_dir)
        
        #establish key parameters, if not supplied by user
        if("character" %in% class(genes_df)){
            genes_df<-
                read_excel_astext(path = genes_df
                                  #,col_types = rep("text",14)
                                  ,sheet = 1)
        }
        
        if(pw_matchup_file=="default"){
            pw_matchup_file<-
                system.file("extdata"
                            ,"pathways_matched_to_sources_current_version.xlsx"
                            ,package="RIGbiopax")
        }
        if(pw_src_file=="default"){
            pw_src_file<-
                system.file("extdata"
                            ,"human_pathways_rbs_curated_current_version.xlsx"
                            ,package="RIGbiopax")
        }
        
        pw_match<-
            read_excel_astext(path = pw_matchup_file
                              #,col_types = rep("text",11)
                              ,sheet = 1) 
        pw_src<-
            read_excel_astext(path = pw_src_file
                              #,col_types = rep("text",4)
                              ,sheet = 1) 
        genes_df<-
            genes_df %>%
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
            genes_df %>%
            RIGbiopax:::internal_analyse_toxdb_biopax_results_df("all")
        
        #summary of dataframe with pathways only from original sources
        summary_list$origsrc<-
            genes_df %>%
            filter(!toxdb.Pathway.ID %in% pw_altsrc) %>%
            RIGbiopax:::internal_analyse_toxdb_biopax_results_df("origsrc")
        
        #summary of dataframe with pathways only from original sources+exact name matches
        summary_list$origsrc_exactmatch<-
            genes_df %>%
            filter(!toxdb.Pathway.ID %in% pw_altsrc) %>%
            filter(pathway.Match.Status=="Exact") %>%
            RIGbiopax:::internal_analyse_toxdb_biopax_results_df("origsrc_exactmatch")
        
        #summary of dataframe with pathways only from original sources+exact name matches+not changed
        summary_list$origsrc_exactmatch_noncur<-
            genes_df %>%
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
            genes_df %>%
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
MAIN_compare_summarize_inxight2biopax<-
    function(wkdir){
        
        #' @title
        #' Compare and Summarize Source BioPAX files to Inxight Pathways List
        #' @description 
        #' Compare genes/proteins extracted from BioPAX to inxight pathways list, and then 
        #' prepare a comprehensive summary such comparison.
        #' @param wkdir Work directory.
        
        #' @author 
        #' Ivan Grishagin
        
        ########### packages
        require(readxl)
        
        ########### workspace with biopax
        biopax_ws<-
            system.file("extdata"
                        ,"biopax_objects_current_version.RData"
                        ,package="RIGbiopax")
        load(biopax_ws
             ,envir = environment())
        
        ########### output dir
        output_dir<-
            file.path(wkdir
                      ,paste0(Sys.Date()
                              ,"_RESULTS_inxight-list_VS_biopax_comparison"))
        
        ########### biopax objects
        biopax_source_names<-
            c("BioCarta"
              ,"KEGG"
              ,"NCI-Nature"
              ,"NetPath"
              ,"Wiki Pathways"
              ,"Science Signaling"
              ,"Reactome"
              ,"RMC"
            )
        biopax_obj<-
            list(biocarta_biopax
                 ,kegg_biopax
                 ,nci_biopax
                 ,netpath_biopax
                 ,wiki_biopax
                 ,scisig_biopax
                 ,reactome_biopax
                 ,rmc_biopax
            )
        
        
        for (ondex in 1:length(biopax_obj)){
            MAIN_compare_toxdb_biopax(work_dir=wkdir
                                      ,source_owl_dir=NULL
                                      ,output_dir=NULL
                                      ,pw_matchup_file="default"
                                      ,toxdb_genes_file="default"
                                      ,source_name=biopax_source_names[ondex]
                                      ,owl_biopax=biopax_obj[[ondex]]
                                      ,verbose=FALSE)
        }
        
        genes_df<-
            merge_gene_comparison_results(dir = output_dir
                                          ,pattern = "genes"
                                          ,filename = NULL)
        pathways_df<-
            merge_gene_comparison_results(dir = output_dir
                                          ,pattern = "pathways"
                                          ,filename = NULL)
        
        #conduct comparisons
        MAIN_biopax_toxdb_comparison_summary(work_dir=wkdir
                                             ,output_dir="./summary"
                                             ,pw_matchup_file="default"
                                             ,pw_src_file="default"
                                             ,genes_df=genes_df
        )
    }

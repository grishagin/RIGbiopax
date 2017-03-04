MAIN_compare_toxdb_biopax <-
    function(work_dir=NULL
             ,source_owl_dir=NULL
             ,output_dir=NULL
             ,pw_matchup_file="default"
             ,toxdb_genes_file="default"
             ,source_name=NULL
             ,owl_biopax=NULL
             ,verbose=FALSE){
        
        #' @title
        #' Compare Proteins Extracted from BioPAX to ToxDB Proteins
        #' @description 
        #' Main function of the package. Extract all proteins ("genes") from a BioPAX source, 
        #' convert their IDs to gene IDs and Symbols, and then compare those to a supplied gene list.
        #' @param work_dir Directory, where all source files/folders are located and will be written to.
        #' @param source_owl_dir Directory with \code{*.owl} files. 
        #' If no directory provided, a prompt will pop up. See \code{load.biopax} function for more information.
        #' @param output_dir .
        #' @param pw_matchup_file File that provides match-up between inxight pathway names and BioPAX pathway names.
        #' @param toxdb_genes_file File that has a table of all genes listed for all pathways. 
        #' @param source_name Name of the BioPAX source.
        #' @param owl_biopax BioPAX object.
        #' @param verbose Logical. Show all relevant progress messages?
        
        #' @author 
        #' Ivan Grishagin
        
        require(RIGessentials)
        
        #establish initial key parameters, if not supplied by user
        if(work_dir=="default"){
            work_dir<-getwd()
        }
        if(pw_matchup_file=="default"){
            pw_matchup_file<-
                system.file("extdata"
                            ,"pathways_matched_to_sources_current_version.xlsx"
                            ,package="RIGbiopax")
        }
        if(toxdb_genes_file=="default"){
            toxdb_genes_file<-
                system.file("extdata"
                            ,"inxight_pathways_current_version.txt"
                            ,package="RIGbiopax")
        }

        prepareSession(work_dir
                       ,nolocale=FALSE)
        
        ########################################################################
        ########################################################################
        #biopax sources
        biopax_source_names<-
            c("BioCarta"          
              ,"KEGG"
              ,"NCI-Nature"
              ,"Reactome"         
              ,"NetPath"         
              ,"Wiki Pathways"
              ,"Science Signaling"
              ,"RMC"
            ) %>%
            sort
        ########################################################################
        ########################################################################
        pkg<-c("plyr"
               ,"dplyr"
               ,"rBiopaxParser"
               ,"mygene"
               ,"readxl"
               ,"openxlsx")
        loadPackages(pkg
                     ,verbose = verbose)
        
        if(is.null(source_name)){
            source_name<-NA
        }
        if(is.null(output_dir)){
            #make directory for output files
            #unless provided
            output_dir<-
                file.path(work_dir
                          ,paste0(Sys.Date()
                                  ,"_RESULTS_inxight-list_VS_biopax_comparison"))
        }
        #create output directory (ignored if exists)
        dir.create(output_dir
                   ,showWarnings = FALSE)
        
        if(!(source_name %in% biopax_source_names)){
            message("Choose BioPAX source name. The choice picker is likely behind your active window.")
            source_name<-
                tkradio_from_vect(biopax_source_names
                                  ,"Select BioPAX Source.")
        }
        
        #load all pathways matched up
        pathways_per_source<-
            read_excel_astext(path = pw_matchup_file
                              #,col_types = rep("text",11)
                              ,sheet = 1) %>%
            dplyr::filter(!is.na(biopax.Pathway.Name)) %>%
            dplyr::select(toxdb.Pathway.ID
                          ,toxdb.Pathway.Name
                          ,biopax.Pathway.Name
                          ,pathway.Match.Status
                          ,Source) %>%
            arrange(Source
                    ,toxdb.Pathway.ID) 
        
        #which source name (biopax) to process
        all_pathways<-
            pathways_per_source %>%
            filter(Source==source_name) 
        all_pathways$toxdb.Pathway.Name<-
            tolower(all_pathways$toxdb.Pathway.Name)
        all_pathways$biopax.Pathway.Name<-
            tolower(all_pathways$biopax.Pathway.Name)
        #choose only pertaining toxdb pathway source and its elements
        toxdb<-
            load_inxight_genes_per_source(source_name=source_name
                                          ,toxdb_genes_file=toxdb_genes_file
                                          ,all_pathways=all_pathways)
        
        #load list of pathways and components from biopax
        if(is.null(owl_biopax)){
            owl_biopax<-
                load.biopax(source_name=source_name
                            ,source_dir=source_owl_dir)
        }
       
        pw_biopax<-
            load.biopax.pathways(owl_biopax=owl_biopax) 
        pw_biopax$biopax.Pathway.Name<-
            pw_biopax$biopax.Pathway.Name %>%
            trimws
        
        #check if pathways that are supposed to be in biopax, actually are there
        mismatch_pathways<-
            #first check if biopax file has all pathways in pathway assignment table 
            all_pathways[!(all_pathways$biopax.Pathway.Name %in% pw_biopax$biopax.Pathway.Name) |
                             #then check if the big table with genes has those pathways as well
                             !(all_pathways$toxdb.Pathway.Name %in% toxdb$toxdb.Pathway.Name),]
        if(nrow(mismatch_pathways)>0){
            message("Some pathways from the list do not match those in biopax or toxdb gene list!")
            print(mismatch_pathways)
            View(mismatch_pathways)
            stop("Aborting.")
        }
        biopax_pw_missing_from_list<-
            pw_biopax[!(pw_biopax$biopax.Pathway.Name %chin% all_pathways$biopax.Pathway.Name),] %>%
            mutate(pathway.Match.Status="Missing"
                   ,toxdb.Pathway.ID=NA
                   ,toxdb.Pathway.Name=NA
                   ,Source=source_name)
        
        ####################### add pathway IDs ################################
        #add biopax pw ids
        all_pathways_pwid<-
            add.biopax.ids(all_pathways=all_pathways
                           ,pw_biopax=pw_biopax)
        
        #add missing pathways from biopax
        all_pathways_pwid<-
            rbind.data.frame(all_pathways_pwid
                             ,biopax_pw_missing_from_list)
        ####################### add pathway IDs ################################
        
        ####################### find pathway components for each ID ################################
        #for each pathway id (row) find the components for that pathway
        #and add them
        df_pw_proteins<-
            adply(.data=all_pathways_pwid
                  ,.margins=1
                  ,.fun=function(x){
                      add_db_ids(owl_biopax=owl_biopax
                                 ,pw_id=x$biopax.Pathway.ID)
                  }) %>%
            mutate(biopax.Gene.Symbol=NA
                   ,ENTREZID=NA) %>%
            data.frame
        
        #add annotations and filter out only entries with db in keytypes
        df_pw_proteins_annot<-
            add.MULT.symbols.entrezids(df_pw_proteins=df_pw_proteins
                                       ,filter_keytypes = TRUE)
        
        #ensure toxdb and annotated biopax df have same columns
        #and fill some of them
        df_pw_proteins_annot<-
            adjust_columns(df_to_adjust = df_pw_proteins_annot
                           ,df_template = toxdb)
        toxdb<-
            adjust_columns(df_to_adjust = toxdb
                           ,df_template = df_pw_proteins_annot)
        
        toxdb$pathway.Match.Status<-
            all_pathways_pwid$pathway.Match.Status[match(toxdb$toxdb.Pathway.ID
                                                         ,all_pathways_pwid$toxdb.Pathway.ID)]
        ####################### find pathway components for each ID ################################
        
        ####################### correlate which genes match ################################
        #prepare data frame(s) with matching genes, genes present in one, but not the other
        #and vice versa
        comparison_results<-
            adply(.data=as.data.frame(all_pathways_pwid)
                  ,.margins=1
                  ,.fun=function(x){
                      compare_pw_components(biopax_prot=df_pw_proteins_annot
                                            ,toxdb_prot=toxdb
                                            ,dfrow=x
                                            ,output="all"
                      )}
            ) %>%
            #remove duplicated biopax pathways
            internal_rm_duplicated_biopax_pw
        #also clean up pathways dataframe accordingly
        all_pathways_pwid<-
            all_pathways_pwid %>%
            filter(biopax.Pathway.ID %in% comparison_results$biopax.Pathway.ID)
        ####################### correlate which genes match ################################
        
        ####################### output #################################
        #output the present and absent pathways into respective files
        openxlsx:::write.xlsx(all_pathways_pwid
                              ,file=
                                  file.path(output_dir
                                            ,paste(Sys.Date()
                                                   ,"pathways"
                                                   ,source_name
                                                   ,"inxight-list_VS_biopax.xlsx"
                                                   ,sep="_"))
                              ,col.names=TRUE
                              ,row.names=FALSE)
        
        #output the comparison results
        openxlsx:::write.xlsx(comparison_results
                              ,file=
                                  file.path(output_dir
                                            ,paste(Sys.Date()
                                                   ,"genes"
                                                   ,source_name
                                                   ,"inxight-list_VS_biopax.xlsx"
                                                   ,sep="_"))
                              ,col.names=TRUE
                              ,row.names=FALSE
                              ,keepNA=TRUE)
        rm(owl_biopax)
        if(verbose){
            invisible(tkmessageBox(message = "All done!"
                                   ,icon = "info"
                                   ,type = "ok")
                      )
        }
                
    }

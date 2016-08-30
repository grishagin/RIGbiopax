MAIN_compare_toxdb_biopax <-
    function(work_dir="D:/Dropbox/Rancho/NCATS/ToxDB"
             ,source_owl_dir=NULL
             ,output_dir=NULL
             ,pw_matchup_file="./_source_files/pathways_matched_to_sources_v013.xlsx"
             ,toxdb_genes_file="./_source_files/toxdb_pathways_15Jun_edIG_2016-08-04.txt"
             ,source_name=NULL){
        
        require(RIGessentials)
        
        
        # work_dir="D:/Dropbox/Rancho/NCATS/ToxDB/"
        # pw_matchup_file="./_source_files/pathways_matched_to_sources_v013.xlsx"
        # toxdb_genes_file="./_source_files/toxdb_pathways_15Jun_edIG_2016-08-04.txt"
        # source_name=NULL
        # source_owl_dir=NULL
        
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
              ,"Ayesha"
              ,"Evgeny"
            )
        ########################################################################
        ########################################################################
        pkg<-c("plyr"
               ,"dplyr"
               ,"rBiopaxParser"
               ,"mygene"
               ,"readxl"
               ,"openxlsx")
        loadPackages(pkg)
        
        if(is.null(source_name)){
            source_name<-NA
        }
        if(is.null(output_dir)){
            #make directory for output files
            #unless provided
            output_dir<-
                file.path(work_dir
                          ,paste0(Sys.Date()
                                  ," RESULTS toxdb-vs-biopax genes comparison"))
            #create said directory
            dir.create(output_dir
                       ,showWarnings = FALSE)
        }
        if(!(source_name %in% biopax_source_names)){
            message("Choose BioPAX source name. The choice picker is likely behind your active window.")
            source_name<-
                tkradio_from_vect(biopax_source_names
                                  ,"Select BioPAX Source.")
        }
        
        #load all pathways matched up
        pathways_per_source<-
            read_excel(path = pw_matchup_file
                       ,col_types = rep("text",11)
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
            load.toxdb.genes.per.source(source_name=source_name
                                        ,source_dir=work_dir
                                        ,toxdb_genes_file=toxdb_genes_file
                                        ,all_pathways=all_pathways)
        
        #load list of pathways and components from biopax
        owl_biopax<-
            load.biopax(source_name=source_name
                        ,source_dir=source_owl_dir)
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
            )
        ####################### correlate which genes match ################################
        
        ####################### output #################################
        #output the present and absent pathways into respective files
        openxlsx:::write.xlsx(all_pathways_pwid
                              ,file=
                                  file.path(output_dir
                                            ,paste(Sys.Date()
                                                   ,"pathways"
                                                   ,source_name
                                                   ,"toxdb_biopax.xlsx"
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
                                                   ,"toxdb_biopax.xlsx"
                                                   ,sep="_"))
                              ,col.names=TRUE
                              ,row.names=FALSE
                              ,keepNA=TRUE)
    }

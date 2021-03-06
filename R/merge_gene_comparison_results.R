merge_gene_comparison_results<-
  function(dir=NULL
           ,pattern=c("genes","pathways")
           ,filename=NULL){
    
    #' @export
    #' @title
    #' Merge Gene Comparison Results
    #' @description 
    #' Merge multiple files that contain gene comparison results by pathway into one file.
    #' @param dir Directory with gene comparison results (in *.xlsx format).
    #' @param pattern Which file name pattern to look for?
    #' @param filename Output filename (optional).
    
    #' @author 
    #' Ivan Grishagin
    
    require(readxl)
    require(dplyr)
    require(openxlsx)
    
    options(stringsAsFactors = FALSE)
    
    files<-list.files(dir,
                      full.names = TRUE,
                      pattern=pattern)
    if(length(files)<1){
      stop("merge_gene_comparison_results: Haven't found any files! Stopping.")
    }
    #ensure pattern is of length 1
    pattern<-
      pattern[1]
    
    #list of dataframes for each excel file
    #which are then merged into one big df
    big_df<-
      lapply(files
             ,read_excel_astext) %>%
      do.call(rbind,.) 
    
    # deprecated
    # if(all(c("ENTREZID","toxdb.Gene.ID","biopax.Gene.ID.Type") %chin% 
    #    colnames(big_df))
    # ){
    #     .genematch<-
    #         rep(NA,nrow(big_df))
    #     .genematch[which(big_df$ENTREZID==big_df$toxdb.Gene.ID)]<-
    #         "1_match"
    #     .genematch[which(is.na(big_df$biopax.Gene.ID.Type) & !is.na(big_df$toxdb.Gene.ID))]<-
    #         "2_toxdb"
    #     .genematch[which(!is.na(big_df$biopax.Gene.ID.Type) & is.na(big_df$toxdb.Gene.ID))]<-
    #         "3_biopax"
    #     .genematch[which(is.na(big_df$biopax.Gene.ID.Type) & is.na(big_df$toxdb.Gene.ID))]<-
    #         "4_genes_not_found"
    #     big_df <-
    #         big_df%>%
    #         mutate(.genematch=.genematch) 
    # }
    
    #look at colnames to arrange by
    arrng_cols<-
      c("Source"
        ,"toxdb.Pathway.ID"
        ,".genematch")
    #see which of them are actually available
    arrng_cols<-
      arrng_cols[arrng_cols %chin% colnames(big_df)]
    
    #arrange the big_df by those cols
    if(length(arrng_cols)>0){
      big_df<-
        big_df %>%
        arrange_(.dots=arrng_cols)
    } else {
      warning("merge_gene_comparison_results: no columns to arrange by, returning unsorted dataframe.")
    }
    
    if(is.null(filename)){
      filename<-
        paste0(Sys.Date()
               ,"_ALL_SOURCES_inxight-list_VS_biopax_"
               ,pattern
               ,".xlsx")
    }
    #write to excel
    openxlsx::write.xlsx(big_df
                         ,file.path(dir
                                    ,filename)
                         ,col.names = TRUE
                         ,row.names = FALSE)
    
    
    return(big_df)
  }
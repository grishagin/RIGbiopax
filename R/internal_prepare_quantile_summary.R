internal_prepare_quantile_summary<-
  function(summary_df
           ,cols=c("matches_toxdb,%"
                   ,"misses_toxdb,%"
                   ,"matches_biopax,%" 
                   ,"extra_biopax,%" )){
    
    #' @keywords internal
    
    if(any(!cols %in% colnames(summary_df))){
      stop("Some of the specified colnames are not among colnames(summary_df)! Aborting.")
    }
    
    #prepare quantile summary
    quantile_summary_df<-
      summary_df[,cols] %>%
      lapply(quantile
             ,seq(0,1,0.1)) %>%
      lapply(round
             ,0) %>%
      do.call(rbind,.) %>%
      cbind.data.frame(Column=rownames(.)
                       ,.)
    
    return(quantile_summary_df)
  }
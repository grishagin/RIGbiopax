adjust_columns <-
function(df_to_adjust=NULL
             ,df_template=NULL){
        if (is.null(df_to_adjust) |
            (is.null(df_template))){
            stop("Not all mandatory parameters provided.")
        }
        #add columns to one df based on the other (template)
        #such that the first would incorporate all from the second
        unique_colnames<-
            unique(c(colnames(df_to_adjust)
                     ,colnames(df_template)))
        
        df_to_adjust[,unique_colnames[!(unique_colnames %in% colnames(df_to_adjust))]]<-
            NA
        return(df_to_adjust)
    }

expand.df.via.split.col <-
    function(dFrame
             ,colsToSplit
             ,patternToSplit=",|;|\\|"){
        #called by add.MULT.symbols.entrezids
        
        if(any(!colsToSplit %in% colnames(dFrame))){
            stop("expand.df.via.split.col: Some colsToSplit values were not found among colnames!")
        }
        
        for(colindex in 1:length(colsToSplit)){
            #which rows are affected in the first column
            rowsAffected<-
                grep(patternToSplit
                     ,dFrame[,colsToSplit[colindex]])
            
            #if any rows are affected, proceed
            if(length(rowsAffected)>0){
                #split the affected rows in the first column by pattern
                splitList<-
                    dFrame[rowsAffected,colsToSplit[colindex]] %>% 
                    strsplit(split=patternToSplit)
                
                #make a vector of how many replicates of each row to take
                rowReplicates<-
                    rep(1,nrow(dFrame))
                rowReplicates[rowsAffected]<-
                    sapply(splitList,length)
                
                finalRowsToTake<-
                    rep(1:nrow(dFrame)
                        ,rowReplicates)
                dFrame<-
                    dFrame[finalRowsToTake,]
                rowsAffected_2<-
                    grep(patternToSplit
                         ,dFrame[,colsToSplit[colindex]])
                
                #replace the values
                dFrame[rowsAffected_2,colsToSplit[colindex]]<-
                    unlist(splitList)
            }
        }
        
        return(dFrame)
    }

expand.df.via.split.col <-
function(dFrame=NULL
             ,colToSplit
             ,patternToSplit=",|;|\\|"){
        #which rows are affected in the first column
        rowsAffected<-
            grep(patternToSplit
                 ,dFrame[,colToSplit[1]])
        
        #if any rows are affected, proceed
        if(length(rowsAffected)>0){
            #split the affected rows in the first column by pattern
            splitList<-
                dFrame[rowsAffected,colToSplit[1]] %>% 
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
                     ,dFrame[,colToSplit[1]])
            
            #replace the values
            dFrame[rowsAffected_2,colToSplit]<-
                unlist(splitList)
        }
        return(dFrame)
    }

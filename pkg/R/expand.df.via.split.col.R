expand.df.via.split.col <-
    function(dFrame
             ,colsToSplit
             ,patternToSplit="\\|"
             ,at_once=TRUE){
        #called by add.MULT.symbols.entrezids
        
        #ensure it's not a tibble
        dFrame<-
            dFrame %>%
            as.data.frame
        
        if(at_once){
            #find affected rows
            rowsAffected<-
                lapply(colsToSplit
                       ,FUN=function(colN){
                           rowN<-
                               grep(patternToSplit
                                    ,dFrame[,colN])
                       }) %>%
                unlist %>%
                unique
            
            #if none affected, return original dataframe
            if (length(rowsAffected)<1){
                return(dFrame)
            }
            
            #for each row do the following
            dFrame_proc<-
                dFrame[rowsAffected,] %>%
                adply(.margin=1
                      ,.fun=function(dfrow){
                          #split every desired column by pattern
                          splitvals_list<-
                              lapply(colsToSplit
                                     ,FUN=function(colN){
                                         unlist(strsplit(dfrow[,colN]
                                                         ,split=patternToSplit))
                                     })
                          #get max of the lengths of every vector in list
                          veclens<-
                              sapply(splitvals_list
                                     ,length)
                          
                          #extend all vectors of length one 
                          #to max vector length (on the list)
                          #ALSO check if length of every vector has the same length
                          splitvals_list<-
                              splitvals_list %>%
                              lapply(FUN=function(vect,N){
                                  if(length(vect)==1){
                                      vect<-
                                          rep(vect,N)
                                  } else if(length(vect)!=N){
                                      message("Problem: max length is\n"
                                              ,N
                                              ,", and vector is "
                                              ,paste(vect,collapse = " "))
                                      print(dfrow)
                                      print(splitvals_list)
                                      stop("expand.df.via.split.col: Split value vectors have unequal lengths!\n")
                                  }
                                  return(vect)
                              }
                              ,N=max(veclens))
                          
                          #now all of the vectors on the list
                          #should have the same lenghts
                          dfrow<-
                              dfrow[rep(1,max(veclens)),]
                          dfrow[,colnames(dfrow) %in% colsToSplit]<-
                              splitvals_list
                          
                          return(dfrow)
                      }) %>%
                #add the unaffected portion of the data frame
                rbind.data.frame(dFrame[-rowsAffected,])
        } else {
            #if it is to be applied iteratively
            #i.e. one column after the other
            dFrame_proc<-
                dFrame
            for(index in 1:length(colsToSplit)){
                dFrame_proc<-
                    expand.df.via.split.col(dFrame_proc
                                            ,colsToSplit=colsToSplit[index]
                                            ,patternToSplit=patternToSplit
                                            ,at_once=TRUE)
            }
        }
        
        return(dFrame_proc)
    }

        
        
        ############################## old version
        # if(any(!colsToSplit %in% colnames(dFrame))){
        #     stop("expand.df.via.split.col: Some colsToSplit values were not found among colnames!")
        # }
        # 
        # for(colindex in 1:length(colsToSplit)){
        #     #which rows are affected in the first column
        #     rowsAffected<-
        #         grep(patternToSplit
        #              ,dFrame[,colsToSplit[colindex]])
        #     
        #     #if any rows are affected, proceed
        #     if(length(rowsAffected)>0){
        #         #split the affected rows in the first column by pattern
        #         splitList<-
        #             dFrame[rowsAffected,colsToSplit[colindex]] %>% 
        #             strsplit(split=patternToSplit)
        #         
        #         #make a vector of how many replicates of each row to take
        #         rowReplicates<-
        #             rep(1,nrow(dFrame))
        #         rowReplicates[rowsAffected]<-
        #             sapply(splitList,length)
        #         
        #         finalRowsToTake<-
        #             rep(1:nrow(dFrame)
        #                 ,rowReplicates)
        #         dFrame<-
        #             dFrame[finalRowsToTake,]
        #         rowsAffected_2<-
        #             grep(patternToSplit
        #                  ,dFrame[,colsToSplit[colindex]])
        #         
        #         #replace the values
        #         dFrame[rowsAffected_2,colsToSplit[colindex]]<-
        #             unlist(splitList)
        #     }
        # }
        # 
        # return(dFrame)

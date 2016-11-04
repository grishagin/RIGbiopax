splitbiopax_from_dt<-
    function(biopax
             ,pwtoextract=NULL
             ,write_to_files=FALSE
             
    ){
        #extract specified (or all) pathways from the biopax along with
        #all of their referenced instances 
        require(RIGbiopax)
        
        #pathways in the biopax
        bppws<-
            listPathways(biopax = biopax) %>%
            .$id
        
        if(is.null(pwtoextract)){
            pwtoextract<-
                bppws
        } else {
            pw_not_present<-
                pwtoextract[!pwtoextract %in% bppws]
            if(length(pw_not_present)>0){
                msg<-
                    paste0("Some pathways could not be found in the biopax table:\n"
                           ,paste(pw_not_present
                                  ,collapse="\n"))
                message(msg)
                #take only those pathways which are present in biopax
                pwtoextract<-
                    pwtoextract[pwtoextract %in% bppws]
            }
        }
        #extract and store all biopax objects as a list
        splitbiopax_list<-
            pwtoextract %>%
            lapply(FUN=function(pwid){
                msg<-
                    paste0("Processing pathway "
                           ,pwid
                           ," ("
                           ,match(pwid,pwtoextract)
                           ," of "
                           ,length(pwtoextract)
                           ,")"
                           )
                message(msg)
                #select all instances referenced by a given pathway
                tempdf<-
                    selectInstances(biopax = biopax
                                    ,id = pwid
                                    ,includeReferencedInstances = TRUE
                                    ,biopaxlevel = 3)
                #if not write to files,
                #filename is NULL
                #so the function will only return a biopax object
                if(!write_to_files){
                    tempfilename<-
                        NULL
                } else {
                    tempfilename<-
                        paste0(pwid
                               ,".owl")
                }
                tempbp<-
                    biopax_from_dt(tempdf
                                   ,filename=tempfilename)
                
                return(tempbp)
            })
        names(splitbiopax_list)<-
            pwtoextract
        
        return(splitbiopax_list)
    }
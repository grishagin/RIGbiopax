splitbiopax_from_dt<-
    function(biopax
             ,pwtoextract=NULL
             ,write_to_files=FALSE
             ,unifyids=TRUE
             
    ){
        #' @title
        #' Convert BioPAX Object to List of BioPAX Objects (by Pathway)
        #' @description 
        #' Splits a BioPAX object into a list of BioPAX objects by pathway. 
        #' Has an option of recording all of these BioPAX objects to files.
        #' @details 
        #' Returns a list of biopax objects.
        #' @param biopax BioPAX object.
        #' @param pwtoextract Vector of pathway ids. If \code{NULL} (default), all pathways are extracted.
        #' @param write_to_files Logical. Write to files? Default is \code{FALSE}.
        #' @param unifyids Logical. Unify IDs by making all elements of one class receive ids in the form "ClassNUMBER", 
        #' with NUMBER ranging from one to number of elements of such class? Default is \code{TRUE}.
.
        
        #' @author 
        #' Ivan Grishagin

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
                if(unifyids){
                    tempdf<-
                        tempdf %>%
                        unify_biopax_ids(idtag=NULL)
                }
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
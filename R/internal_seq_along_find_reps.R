internal_seq_along_find_reps<-
  function(vect){
    
    #' @keywords internal
    
    #returns a vector where all unique positions are 1's
    #and their duplicates are 2's, and so on
    #a,b,a,a,c,b,d -converted_into-> 1,1,2,3,1,2,1
    toreturn<-rep(NA,length(vect))
    for(unid in unique(vect)){
      toreturn[vect %in% unid]<-
        1:sum(vect %in% unid)
    }
    return(toreturn)
  }
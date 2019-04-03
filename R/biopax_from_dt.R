biopax_from_dt<-
  function(dTable
           ,filename=NULL
           ,encoding="UTF-8"){
    #' @export
    #' @title
    #' Prepare BioPAX Object from BioPAX-Style DataTable
    #' @description 
    #' Prepares a BioPAX object from a BioPAX-style datatable.
    #' @details 
    #'
    #' @param dTable BioPAX-style datatable.
    #' @param filename (optional) Output filename. 
    #' If provided, BioPAX object will also be written to file.
    #' @param encoding (currently not used).
    
    #' @author 
    #' Ivan Grishagin
    
    #makes biopax object from a biopax-style data table
    #and optionally writes it to file
    
    #create new_biopax
    new_biopax<-
      createBiopax()
    #add data table
    new_biopax$dt<-
      rbind(new_biopax$dt
            ,dTable)
    class(new_biopax$dt)<-
      c("biopax_df"
        ,class(new_biopax$dt))
    if(!is.null(filename)){
      #write to file
      writeBiopax_Rancho(biopax = new_biopax
                         ,filename = filename
                         ,overwrite = TRUE
                         ,biopaxlevel =3)
    }
    return(new_biopax)
  }
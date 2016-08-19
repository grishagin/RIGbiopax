biopax_from_dt<-
    function(dTable
             ,filename=NULL
             ,encoding="ISO-8859-1"){

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
                               ,file = filename
                               ,overwrite = TRUE
                               ,biopaxlevel =3
                               ,encoding = encoding)
        }
        return(new_biopax)
    }
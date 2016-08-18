.checkValidity_Rancho<-
    function (biopax) 
    {
        #copy from rBiopaxParser
        
        if (!any(grepl("biopax", class(biopax)))) 
            stop("Supplied biopax object doesnt seem to be of class biopax!")
        if (nrow(biopax$dt) < 1) 
            stop("Internal data.frame of supplied biopax object seems to be empty!")
        if (ncol(biopax$dt) != 6) 
            stop("Internal data.frame of supplied biopax object seems to be invalid!")
    }
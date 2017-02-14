update_RIGbiopax <-
    function () 
    {
        #' @title
        #' Update RIGbiopax from GitHub
        #' @description 
        #' Updates the package. It may be necessary to restart the R session after the update.
        
        #' @author 
        #' Ivan Grishagin
        
        
        #unloadNamespace(ns = "RIGbiopax")
        devtools::unload(pkg = inst("RIGbiopax"))
        devtools::install_github("grishagin/RIGbiopax", subdir = "pkg")
    }
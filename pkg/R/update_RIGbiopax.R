update_RIGbiopax <-
    function () 
    {
        unloadNamespace(ns = "RIGbiopax")
        remove.packages("RIGbiopax")
        devtools::install_github("grishagin/RIGbiopax", subdir = "pkg")
    }
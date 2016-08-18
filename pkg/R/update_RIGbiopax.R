update_RIGbiopax <-
    function () 
    {
        unloadNamespace(ns = "RIGbiopax")
        devtools::install_github("grishagin/RIGbiopax", subdir = "pkg")
    }
striphash<-
  function (vect) {
    #' @export
    #' @title
    #' Remove # from String
    #' @description
    #' Removes the first \code{#} from each element of a string vector.
    #' @param vect String vector.
    #'
    #' @author
    #' Ivan Grishagin
    
    sub("#"
        , ""
        , vect)
  }
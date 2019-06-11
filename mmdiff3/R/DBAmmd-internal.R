# Internal function to set distances to DBAmmd Object
#
# This function is used internally in the S4 method DBAmmd, and should not be
# called by the user in general.
#

setDists <- function(x,DISTs,name) {
  D <- Dists(x)
  if (!is.null(D[[name]])){
    warning(paste('DISTs (', name, ' already exists.)') )
    D[[name]] <- cbind(D[[name]],DISTs)
    x@DISTs <- D
  } else{
    D[[name]] <- DISTs
    x@DISTs <- D
  }
  return(x)}

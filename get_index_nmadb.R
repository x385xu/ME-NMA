library(netmeta)
library(nmadb)
library(dplyr)


# Extract index of datasets in nmadb based on the measure
#  measure: string, like "risk ratio" or "odds ratio"

get_index_nmadb <- function(dat = dat_nmadb,
                            measure) {
  # Extract index with the according measure
  ind_measure <- which(dat_nmadb$Effect.Measure == measure)
  
  twoarm <- c()
  multiarm <- c()

  for (i in ind_measure) {
    # return null if runnetmeta cannot get the dataset
    net <- tryCatch({runnetmeta(dat_nmadb$recid[i])}, 
                    error = function(e) {return(NULL)})
    
    # Skip if net is not a list
    if (!is.list(net)) next
    
    has_multi <- any(net$multiarm) # is FALSE if every entry is FALSE (no multiarm studies)
    
    if (has_multi) {
      multiarm <- c(multiarm, i)
    } else {
      twoarm <- c(twoarm, i)
    }
  }
  n <- max(length(twoarm), length(multiarm))
  twoarm_p   <- c(twoarm,   rep(NA_integer_, n - length(twoarm)))
  multiarm_p <- c(multiarm, rep(NA_integer_, n - length(multiarm)))
  
  return(data.frame(twoarm = twoarm_p, multiarm = multiarm_p))
}



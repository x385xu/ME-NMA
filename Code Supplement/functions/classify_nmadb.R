#' @description
#' Classify networks in nmadb by effect measure and study design
#'
#' This function iterates over a network meta-analysis database and classifies
#' each network according to its effect measure and whether it consists solely
#' of two-arm studies or includes multi-arm trials.
#'
#' @param dat Data frame containing the network meta-analysis database.
#'   Must include:
#'   \itemize{
#'     \item \code{recid}: Unique identifier for each network (used by \code{runnetmeta()}).
#'     \item \code{Effect.Measure}: Effect measure for each network (e.g., "odds ratio",
#'       "risk ratio", "mean difference").
#'   }
#'   Defaults to \code{dat_nmadb}.
#'
#' @return A data frame with one row per successfully processed network, containing:
#' \itemize{
#'   \item \code{index}: Row index of the network in \code{dat}.
#'   \item \code{recid}: Record ID of the network.
#'   \item \code{effect_measure}: Effect measure associated with the network.
#'   \item \code{twoarm_only}: Logical indicator; \code{TRUE} if the network
#'     contains only two-arm studies, \code{FALSE} if at least one multi-arm
#'     study is present.
#' }
#'


classify_nmadb <- function(dat = dat_nmadb) {
  
  out <- data.frame(
    index = integer(),
    recid = character(),
    effect_measure = character(),
    twoarm_only = logical()
  )
  
  for (i in seq_len(nrow(dat))) {
    recid_i <- dat_nmadb$recid[i]
    dat_i <- readByID(recid_i)
    
    net <- tryCatch(
      fit_netmeta(dat_i),
      error = function(e) NULL
    )
    
    if (!is.list(net)) next
    
    has_multi <- any(net$multiarm)
    
    out <- rbind(
      out,
      data.frame(
        index = i,
        recid = dat$recid[i],
        effect_measure = dat$Effect.Measure[i],
        twoarm_only = !has_multi
      )
    )
  }
  
  return(out)
}
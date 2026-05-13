#' @description
#' Fit a network meta-analysis model using the \code{netmeta} package.
#'
#' Parts of this function adapt functionality and workflow ideas from the
#' \code{nmadb} package (GPL-3 licensed), originally developed for fitting
#' network meta-analysis models from the NMA database. This implementation
#' was rewritten to operate on locally saved datasets and to avoid reliance
#' on the archived and no-longer-maintained \code{nmadb} package API.
#'
#' The function supports binary, continuous, rate, and inverse-variance
#' formatted datasets. Arm-based datasets are first converted to
#' pairwise comparisons using \code{meta::pairwise()}, after which the
#' model is fitted using \code{netmeta::netmeta()}.
#'
#' @param indata A single NMA dataset object stored locally. The object 
#'   must contain metadata fields such as \code{type}, \code{format},
#'   \code{effect}, and \code{data}.
#' @param model Character string specifying the model type.
#'   Either \code{"common"} for a common-effect model or
#'   \code{"random"} for a random-effects model.
#'   Defaults to \code{"random"}.
#' @param measure Optional summary measure passed to
#'   \code{netmeta}. If \code{NULL}, the original effect measure
#'   stored in the dataset is used.
#' @param method.tau Character string specifying the estimator used
#'   for the between-study heterogeneity standard deviation
#'   \eqn{\tau}. Passed to \code{netmeta()}.
#'   Defaults to \code{"DL"}.
#' @param details.chkmultiarm Logical indicating whether detailed
#'   diagnostics for multi-arm studies should be returned by
#'   \code{netmeta()}. Defaults to \code{TRUE}.
#' @param ... Additional arguments passed directly to
#'   \code{netmeta::netmeta()}.
#'
#' @return
#' An object of class \code{netmeta}.
#'
#' @seealso
#' \code{\link[netmeta]{netmeta}},
#' \code{\link[meta]{pairwise}}
#'

fit_netmeta <- function(indata,
                        model = "random",
                        measure = NULL,
                        method.tau = "DL",
                        details.chkmultiarm = TRUE,
                        control = NULL,
                        ...) {
  
  if (is.null(indata)) stop("`indata` is NULL.")
  
  type <- indata$type
  format <- as.character(indata$format)
  
  # Some datasets in nmadb are labelled as "other"
  # and are not compatible with netmeta
  if (is.null(measure)) {
    if (indata$effect == "other") {
      stop("Cannot analyze atypical effect measure: `other`.")
    }
    sm <- indata$effect
  } else {
    sm <- measure
  }
  
  # Construct internal data type label
  # e.g. "long_binary", "long_continuous"
  data_type <- if (format != "iv") {
    paste0("long_", type)
  } else {
    "iv"
  }
  
  D <- indata$data
  
  # Wide-format datasets are not currently supported
  if (format == "wide") {
    stop(
      "This function currently expects `long` or `iv` formatted data. ",
      "Convert wide-format datasets before fitting."
    )
  }
  
  common_fit <- model %in% c("fixed", "common")
  random_fit <- model %in% c("random")
  
  # ------------------------------------------------------------
  # Convert arm-based data to pairwise comparisons
  # using meta::pairwise()
  # ------------------------------------------------------------
  if (data_type == "long_binary") {
    
    # Binary outcome data
    pw <- meta::pairwise(
      treat = D$t,
      event = D$r,
      n = D$n,
      studlab = D$id,
      data = D,
      sm = sm,
      allstudies = TRUE
    )
    
  } else if (data_type == "long_continuous") {
    
    # Continuous outcome data
    pw <- meta::pairwise(
      treat = D$t,
      mean = D$y,
      sd = D$sd,
      n = D$n,
      studlab = D$id,
      data = D,
      sm = sm,
      allstudies = TRUE
    )
    
  } else if (data_type == "long_rate") {
    
    # Rate outcome data
    pw <- meta::pairwise(
      treat = D$t,
      event = D$r,
      time = D$time,
      studlab = D$id,
      data = D,
      sm = sm,
      allstudies = TRUE
    )
    
  } else if (data_type == "iv") {

    # Inverse-variance format already contains treatment effects
    # and standard errors, so no pairwise conversion is needed
    fit <- netmeta::netmeta(
      TE = D$effect,
      seTE = D$se,
      treat1 = D$t1,
      treat2 = D$t2,
      studlab = D$id,
      data = D,
      sm = sm,
      common = common_fit,
      random = random_fit,
      method.tau = method.tau,
      details.chkmultiarm = details.chkmultiarm,
      control = control,
      ...
    )
    
    return(fit)
    
  } else {
    stop("Unsupported data type: ", data_type)
  }
  
  fit <- netmeta::netmeta(
    TE = pw$TE,
    seTE = pw$seTE,
    treat1 = pw$treat1,
    treat2 = pw$treat2,
    studlab = pw$studlab,
    data = pw,
    sm = sm,
    common = common_fit,
    random = random_fit,
    method.tau = method.tau,
    details.chkmultiarm = details.chkmultiarm,
    control = control,
    ...
  )
  
  return(fit)
}


#' Simulates an RxODE model
#'
#' @param model the RxODE model
#' @param event_table the RxODE event table
#' @param inits a vector of initial values of the state variables (e.g., amounts in each compartment),
#' and the order in this vector must be the same as the state variables (e.g., PK/PD compartments)
#' @param THETA model THETAS
#' @param sim_times simulation times
#' @param eta model ETAS or NULL if only typical values should be simulated
#' @param iCov dataframe of individual time-independent covariates or NULL if not used
#' @param ... additonal parameters for RxODE 'solve' function
#'
#' @return RxODE simulation result
#' @export
#'
simulate_model <- function(model,
                           event_table,
                           inits,
                           THETA,
                           sim_times,
                           ETA = NULL,
                           iCov = NULL,
                           ...) {

  # get model parameter names that are not THETAS
  param.names <- model$params
  # eleminate thetas
  param.names <- param.names[-which(.grepl_theta(param.names))]

  # replace missing etas (shortcut for simulations only)
  eta.entries <- which(.grepl_eta(param.names))
  n.eta <- length(eta.entries)
  if (n.eta > length(ETA)) {
    n.eta.missing <- n.eta - length(ETA)
    ETA <- c(ETA, rep(0, n.eta.missing))
    message(paste("Added", n.eta.missing, "ETA entries with 0"))
  }

  # eleminate etas from parameter names
  if (n.eta > 0) {
    param.names <- param.names[-eta.entries]
  }

  # check for time-independent covariates
  if (!is.null(iCov)) {
    i.cov.names <- names(iCov)
    if (length(param.names) == 0)
      stop("iCov defined but no covariate found in model", call. = FALSE)

    icov.idx <- which(i.cov.names %in% param.names)
    if (length(icov.idx) != length(i.cov.names))
      stop("Could not match iCov names with model covariates", call. = FALSE)

    param.names <- param.names[-icov.idx]
  }

  ev_t <- event_table$copy()
  ev_t$add.sampling(sim_times)

  # solve
  result <- model %>% RxODE::rxSolve(params = iCov,
                                     events = ev_t,
                                     inits = inits,
                                     theta = THETA,
                                     eta = ETA,
                                      ...)
  return(result)
}



# map objective function
# ETA -> Eta of PK-Parameter of the individuum (row vector)
# I_OMEGA -> INVERSE of IIV of the PK-Parameter (VARIANCE)
# C_i -> Plasma concentrations of the indiviuum
# C_m -> Plasma concentraions of the model prediction
# SIGMA -> Residual variance of the model (VARIANCE)
# !! Only works with proporional error model !!


#' Objective function of the empicical bayes estimate (EBE)
#'
#' @param ETA ETAs of PK/PD-Parameter of the individual as row vector
#' @param I_OMEGA Inverve of the IIV of the PK-Parameter (interindividual variance)
#' @param C_i observed values of the individual
#' @param C_m simulated values for the individual
#' @param SIGMA Residual variance of the model
#' @param error_model The residual error model/function
#'
#' @return Objective function value
#'
map_obj_fn <- function(ETA, I_OMEGA, C_i, C_m, SIGMA,
                       error_model = c("proportional", "additive")) {

  error_model <- match.arg(error_model)

  # see e.g
  #   https://synapse.koreamed.org/DOIx.php?id=10.4196/kjpp.2012.16.2.97
  #   Eq 2

  # ETAs
  left_side = 0
  if (length(ETA) > 0)
    left_side = t(ETA) %*% I_OMEGA %*% ETA


  # DVs
  right_side = 0
  if (length(C_i) > 0) {
    if (error_model == "proportional") {

      # parital div of the error function Y = F + F * EPS(i)
      # dY/dEPS(i) = F
      # -> Cov-Matrix
      # prop.error = SIGMA
      # sigma_sq = diag(Cov)**2 * prop.error
      s_sq <- SIGMA * (C_m**2)
    } else if (error_model == "additive") {

      # Y = F + EPS(i)
      # dY/dEPS(i) = 1
      s_sq <- SIGMA
    } else {
      stop("error_model not implemented")
    }

    res = (C_i - C_m)**2 / s_sq
    right_side = sum(res + log(s_sq))
  }

  res <- left_side + right_side
  return(res)
}



#' Creates an MAP objective function
#'
#' @param individual_df NONMEM dataframe with individual measurements
#' @param model the RxODE model
#' @param event_table the RxODE event table
#' @param inits a vector of initial values of the state variables (e.g., amounts in each compartment),
#' and the order in this vector must be the same as the state variables (e.g., PK/PD compartments)#'
#' @param THETA Typical values THETA
#' @param I_OMEGA Inverse of model OMEGA as matrix
#' @param SIGMA residual variance SIGMA of the model
#' @param iCov dataframe of individual time-independent covariates or NULL if not used
#' @param error_model The residual error model/function
#' @param ... additonal parameters for RxODE 'solve' function
#'
#' @return the MAP objective function
#' @export
#'
create_obj_fn <- function(individual_df,
                          model,
                          event_table,
                          inits,
                          THETA,
                          I_OMEGA,
                          SIGMA,
                          iCov = NULL,
                          error_model = c("proportional", "additive"),
                          ...) {

  # proportial: Y = F * (1 + e)
  # additive  : Y = F + e
  error_model <- match.arg(error_model)

  pred.col =  "IPRED"
  if (!(pred.col %in% model$lhs))
    stop("Could not find lhs variable IPRED")

  DVs <- obs_values(individual_df, include_time = FALSE)
  obs_times <- obs_times(individual_df, include_events = FALSE)

  obj_fn <- function(eta) {

    simulation <- simulate_model(model,
                                 event_table,
                                 inits,
                                 THETA,
                                 obs_times,
                                 ETA = eta,
                                 iCov = iCov, ...)

    IPRED <- simulation[["IPRED"]]
    res <- map_obj_fn(eta, I_OMEGA, DVs, IPRED, SIGMA, error_model)
    return(res)
  }

  return(obj_fn)
}


#' MAP estimation for model ETAs for an individual
#'
#' @param individual_df NONMEM dataframe with measured data of the individual for the estimation
#' @param model RxODE model
#' @param event_table RxODE event table
#' @param inits a vector of initial values of the state variables (e.g., amounts in each compartment),
#' and the order in this vector must be the same as the state variables (e.g., PK/PD compartments)#'
#' @param THETA values of structural THETA parameters of the model.
#' @param OMEGA variance of intersubject errors. Must be a matrix or vector.
#' @param SIGMA variance
#' @param iCov dataframe of individual time-independent covariates or NULL if not used
#' @param error_model The residual error model/function
#' @param fit_options Options for the optim fitting routine
#' @param ... additonal parameters for RxODE 'solve' function
#'
#' @return List with optimized parameters, hessian, information about convergence and standard errors of the fit.
#' @export
#'
map_estimate <- function(individual_df,
                         model,
                         event_table,
                         inits,
                         THETA,
                         OMEGA,
                         SIGMA,
                         iCov = NULL,
                         error_model = c("proportional", "additive"),
                         fit_options = list(method = "L-BFGS-B"),
                         ...) {

  # proportial: Y = F * (1 + e)
  # additive  : Y = F + e
  error_model <- match.arg(error_model)


  # calculate the inverse of OMEGA
  if (is.null(dim(OMEGA)))
    I_OMEGA <- diag(1/OMEGA)     # vector
  else
    I_OMEGA <- solve(OMEGA)      # matrix

  # create objective function
  obj_fn <- create_obj_fn(individual_df,
                          model,
                          event_table,
                          inits,
                          THETA,
                          I_OMEGA,
                          SIGMA,
                          iCov = iCov,
                          error_model = error_model,
                          ...)

  # fit
  start_params <-  rep(0, nrow(I_OMEGA))
  fit <- do.call(optim, c(list(par = start_params,
                               fn = obj_fn,
                               hessian = TRUE),
                          fit_options))

  if (fit$convergence != 0)
    warning(paste("Fit did not converge:", fit$message))

  return(list(par = fit$par,
              hess = fit$hessian,
              converged = (fit$convergence == 0),
              se = sqrt(diag(2 * abs(solve(fit$hessian))))))

}

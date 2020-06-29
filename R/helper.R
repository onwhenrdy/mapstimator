#' @importFrom magrittr "%>%"


# grepl ETAS from vector of names
.grepl_eta <- function(names) {
  grepl("^ETA\\[[0-9]*\\]$", names)
}

# grepl THETAS from vector of names
.grepl_theta <- function(names) {
  grepl("^THETA\\[[0-9]*\\]$", names)
}

# removes event rows from a NONMEM data.frame
.remove_events <- function(df) {

  event_lines <- which(df[["EVID"]] != 0)
  if (length(event_lines) > 0)
    df <- df[-event_lines,]

  return(df)
}


#' Splits dataframe into list of dataframes by id column 'ID'
#'
#' @param df THe data.frame
#'
#' @return A named list of data.frames with names = ID
#' @export
#'
split_by_id <- function(df) {

  id_col <- "ID"
  ids <- unique(df[[id_col]])
  if (length(ids) < 1 || is.na(ids)) {
    stop("Did not find ids")
  }

  results <- list()
  for (id in ids) {
    idx <- which(df[[id_col]] == id)
    tmp <- list(temp = df[idx,])
    names(tmp) <- id
    results <- append(results, tmp)
  }
  return(results)
}


#' Creates a init vector for RxODE from the model compartment list
#'
#' @param model the RxODE model
#'
#' @return Init vector initialized with 0s
#' @export
#'
create_init <- function(model) {

  comp.names <- model$state
  inits <- rep(0, length(comp.names))
  names(inits) <- comp.names

  return(inits)
}


#' Adds dosing events for a model from a NONMEM dataframe
#' Works for bolus doses and infusions.
#' Function does not modify the input event table.
#'
#' @param event_table The event table that should added to
#' @param model The target model
#' @param df The NONMEM dataframe
#'
#' @return A modified event table with added dosing events
#' @export
#'
add_dosing_events <- function(event_table, model, df) {

  ev_t <- event_table$copy()
  comp.size <- length(model$state)

  event_lines <- which(df[["EVID"]] != 0)
  if (length(event_lines) > 0) {

    event_df <- df[event_lines,]
    for (i in 1:nrow(event_df)) {
      comp.idx <- event_df[["CMT"]][i]

      if (is.na(comp.idx) || comp.idx > comp.size) {
        stop(paste("Compartment", comp.idx, "cannot be found in the model"))
      }

      dose <- event_df[["AMT"]][i]
      if (is.null(dose) || is.na(dose)) {
        stop("Dose is NA or NULL")
      }

      time <- event_df[["TIME"]][i]
      if (is.null(time) || is.na(time)) {
        stop("Time is NA")
      }

      rate <- event_df[["RATE"]][i]

      ev_t$add.dosing(dose = dose, start.time = time, rate = rate)
    }
  }

  return(ev_t)
}


#' Extracts observed times from a NONMEM dataframe
#'
#' @param df the NONMEM dataframe
#' @param include_events TRUE if event times should be included
#'
#' @return A vector of times with unique times. Result will be sorted.
#' @export
#'
obs_times <- function(df, include_events = FALSE) {

  if (!include_events) {
    df <- .remove_events(df)
  }

  result <- df[["TIME"]]
  result <- unique(result)
  result <- sort(result)

  return(result)
}


#' Extracts observed values from a NONMEM dataframe
#'
#' @param df the NONMEM dataframe
#' @param include_time TRUE, if the time column should be extracted as well, else FALSE
#'
#' @return A dataframe with observed values
#' @export
#'
obs_values <- function(df, include_time = FALSE) {

  idx <- c("DV")
  if (include_time) {
    idx <- c("TIME", idx)
  }

  event_lines <- which(df[["EVID"]] != 0)
  df <- df[-event_lines,]

  return(df[idx])
}

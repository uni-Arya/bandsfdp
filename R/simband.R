#' Calculate simultaneous FDP bounds given a set of ordered labels resulting
#' from TDC.
#'
#' In TDC, each hypothesis is attached to a winning score and a label (1 for a
#' target win and -1 for a decoy win). \code{simband()} assumes that the
#' hypotheses are ordered in decreasing order of winning scores (breaking
#' ties at random). In particular, the argument \code{labels} is ordered
#' according to this rule.
#'
#' @param labels A vector of (ordered) labels resulting from TDC.
#' @param gamma The confidence level of the band. Typical values include
#' \code{gamma = 0.05} or \code{gamma = 0.01}.
#' @param type A character string specifying which band to use. Must be one of
#' "stband" or "uniband".
#' @param indices A vector of indices of (ordered) hypotheses associated to TDC.
#' Defaults to \code{1:n}, where \code{n} is the number of hypotheses.
#' @param d_max An optional positive integer specifying the maximum number
#' of decoy wins considered in calculating the bands.
#' @param max_fdp A number specifying the maximum FDP considered by the user in
#' calculating the bands. Used to compute dmax if dmax is set to NULL.
#' @param c Determines the ranks of the target score that are considered
#' winning. Defaults to \code{c = 0.5} for TDC.
#' @param lambda Determines the ranks of the target score that are
#' considered losing. Defaults to \code{lambda = 0.5} for TDC.
#' @param interpolate A boolean indicating whether the bands should be
#' interpolated. Offers a slight boost in performance at the cost of computing
#' power. Defaults to \code{TRUE}.
#'
#' @return A vector of upper prediction bounds on the FDP for each of the
#' specified indices.
#' @export
#'
#' @examples
#' if (requireNamespace("fdpbandsdata", quietly = TRUE)) {
#'   set.seed(123)
#'   labels <- c(
#'     rep(1, 250),
#'     sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.9, 0.1)),
#'     sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.5, 0.5)),
#'     sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.1, 0.9))
#'   )
#'   gamma <- 0.05
#'   simband(labels, gamma, "stband", indices = c(100, 250, 500, 1000))
#' }
simband <- function(labels, gamma, type,
                    indices = 1:length(labels),
                    d_max = NULL, max_fdp = 0.5,
                    c = 0.5, lambda = 0.5,
                    interpolate = TRUE) {
  if (any(indices < 0)) {
    stop("Negative threshold detected.")
  }
  if (!type %in% c("stband", "uniband")) {
    stop("Invalid type argument in simband().")
  }
  if (requireNamespace("fdpbandsdata", quietly = TRUE)) {
    n <- length(labels)
    B <- c / (1 - lambda)

    # Name of the table to search for in fdpbandsdata
    name_search <- function(num) {
      format(ceiling(num * 10^5 - 1e-12) / 10^5, nsmall = 5)
    }

    if (type == "stband") {
      # Find the table which contains z (to be used in xi_sb())
      table_name <- paste0(
        "qtable_c",
        name_search(c), "_lam", name_search(lambda)
      )
      table <- get(table_name, envir = getNamespace("fdpbandsdata"))

      xi_sb <- function(d, d_max) {
        if (d + 1 > d_max) {
          return(Inf)
        }
        max(0, floor(z * sqrt(B * (1 + B) * (d + 1)) + B * (d + 1) + 1e-9))
      }
    } else {
      # Find the table which contains u (to be used in xi_ub())
      table_name <- paste0(
        "utable_c",
        name_search(c), "_lam", name_search(lambda)
      )
      table <- get(table_name, envir = getNamespace("fdpbandsdata"))

      xi_ub <- function(d, d_max) {
        if (d + 1 > d_max) {
          return(Inf)
        }
        floor(stats::qnbinom(1 - u, d + 1, prob = (1 / (1 + B))) + 1e-9)
      }
    }

    if (!is.null(d_max)) {
      if (d_max <= 0) {
        stop("Argument d_max is non-positive.")
      }
      if (d_max > 5e4) {
        warning("d_max is greater than 50000. Setting d_max = 50000.")
        d_max <- 5e4
      }

      if (type == "stband") {
        z <- table[d_max, paste0(1 - gamma)]
      } else {
        rho <- table[[1]][d_max, paste0(1 - gamma)]
        sigma <- table[[2]][d_max, paste0(1 - gamma)]
        flip <- table[[3]][d_max, paste0(1 - gamma)]
        u <- ifelse(stats::runif(1) < flip, rho, sigma)
      }
    } else {
      d_max <- 0
      res <- 0
      if (type == "stband") {
        while (res < max_fdp) {
          if (d_max + 1 > n) {
            d_max <- n
            break
          }

          if (d_max + 1 > 5e4) {
            warning("Computed d_max is at least 50000. Setting d_max = 50000.")
            d_max <- 5e4
            break
          }

          z <- table[d_max + 1, paste0(1 - gamma)]

          res <- xi_sb(d_max, d_max + 1) / (n - d_max)
          if (res > max_fdp) {
            break
          }
          d_max <- d_max + 1
        }
      } else {
        while (res < max_fdp) {
          if (d_max + 1 > n) {
            d_max <- n
            break
          }

          if (d_max + 1 > 5e4) {
            warning("Computed d_max is at least 50000. Setting d_max = 50000.")
            d_max <- 5e4
            break
          }

          rho <- table[[1]][d_max + 1, paste0(1 - gamma)]
          sigma <- table[[2]][d_max + 1, paste0(1 - gamma)]
          flip <- table[[3]][d_max + 1, paste0(1 - gamma)]
          u <- ifelse(stats::runif(1) < flip, rho, sigma)

          res <- xi_ub(d_max, d_max + 1) / (n - d_max)
          if (res > max_fdp) {
            break
          }
          d_max <- d_max + 1
        }
      }

      if (d_max == 0) {
        stop("Calculated d_max is 0. The bounds over the indices will be 1.")
      }

      if (type == "stband") {
        z <- table[d_max, paste0(1 - gamma)]
      } else {
        rho <- table[[1]][d_max, paste0(1 - gamma)]
        sigma <- table[[2]][d_max, paste0(1 - gamma)]
        flip <- table[[3]][d_max, paste0(1 - gamma)]
        u <- ifelse(stats::runif(1) < flip, rho, sigma)
      }
    }

    decoy_fun <- function(threshold) {
      if (threshold > 0) {
        sum(labels[1:threshold] == -1)
      } else {
        0
      }
    }

    target_fun <- function(threshold) {
      if (threshold > 0) {
        sum(labels[1:threshold] == 1)
      } else {
        0
      }
    }

    if (interpolate) {
      dbar <- function(threshold, type) {
        check_these <- c(
          1,
          setdiff(which(labels[1:threshold] == -1) - 1, 0),
          threshold
        )

        max_value <- 0

        if (type == "stband") {
          for (t in check_these) {
            new_value <- ceiling(
              target_fun(t) - xi_sb(decoy_fun(t), d_max) - 1e-12
            )
            if (new_value > max_value) {
              max_value <- new_value
            }
          }
        } else {
          for (t in check_these) {
            new_value <- ceiling(
              target_fun(t) - xi_ub(decoy_fun(t), d_max) - 1e-12
            )
            if (new_value > max_value) {
              max_value <- new_value
            }
          }
        }

        return(max_value)
      }

      with_interpolation <- function(threshold, type) {
        if (threshold == 0) {
          0
        } else if (decoy_fun(threshold) >= 5e4) {
          1
        } else {
          (target_fun(threshold) - dbar(threshold, type)) /
            max(1, target_fun(threshold))
        }
      }

      return(mapply(with_interpolation, indices, rep(type, length(indices))))
    } else {
      without_interpolation <- function(d, t) {
        if (t == 0) {
          return(0)
        }
        if (type == "stband") {
          min(1, xi_sb(d, d_max) / max(t, 1))
        } else {
          min(1, xi_ub(d, d_max) / max(t, 1))
        }
      }

      decoys <- vapply(indices, function(t) decoy_fun(t), numeric(1))
      targets <- vapply(indices, function(t) target_fun(t), numeric(1))

      return(mapply(without_interpolation, decoys, targets))
    }
  } else {
    stop(
      "Simulataneous bands require precomputed data tables. You may choose",
      " to run devtools::install_github(\"uni-Arya/fdpbandsdata\") to install",
      " these tables."
    )
  }
}

#' Calculate the standardized band on TDC's list of discoveries
#'
#' In TDC, each hypothesis is attached to a winning score and a label (1 for a
#' target win and -1 for a decoy win). \code{stband()} assumes that the
#' hypotheses are ordered in decreasing order of winning scores (breaking
#' ties at random). In particular, the argument \code{labels} is ordered
#' according to this rule.
#'
#' @param thresholds A vector of rejection thresholds resulting from TDC.
#' @param labels A vector of (ordered) labels resulting from TDC.
#' @param alpha The FDR threshold.
#' @param gamma The confidence level of the band. Typical values include
#' \code{gamma = 0.05} or \code{gamma = 0.01}.
#' @param c Determines the ranks of the target score that are considered
#' winning. Defaults to \code{c = 0.5} for TDC.
#' @param lambda Determines the ranks of the target score that are
#' considered losing. Defaults to \code{lambda = 0.5} for TDC.
#' @param n The number of hypotheses. Defaults to the length of \code{labels}.
#' @param interpolate A boolean indicating whether the bands should be
#' interpolated. Offers a slight boost in performance at the cost of computing
#' power. Defaults to \code{TRUE}.
#'
#' @return A vector of upper prediction bounds on the FDP for each rejection
#' threshold.
#' @export
#'
#' @examples
#' if (requireNamespace("fdpbandsdata", quietly = TRUE)) {
#'   set.seed(123)
#'   thresholds <- c(250, 500, 750, 1000)
#'   labels <- c(
#'     rep(1, 250),
#'     sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.9, 0.1)),
#'     sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.5, 0.5)),
#'     sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.1, 0.9))
#'   )
#'   alpha <- 0.05
#'   gamma <- 0.05
#'   stband(thresholds, labels, alpha, gamma)
#' }
stband <- function(thresholds, labels,
                   alpha, gamma,
                   c = 0.5, lambda = 0.5,
                   n = length(labels),
                   interpolate = TRUE) {
  if (any(thresholds < 0)) {
    stop("Negative threshold detected.")
  }
  if (requireNamespace("fdpbandsdata", quietly = TRUE)) {
    B <- c / (1 - lambda)
    d_max <- max(1, floor(alpha * (n + 1) / (alpha + B) + 1e-12))

    if (d_max > 5e4) {
      warning("d_max is greater than 50000. Setting d_max = 50000.")
      d_max <- 5e4
    }

    # Look up the appropriate quantile table
    name_search <- function(num) {
      format(ceiling(num * 10^5 - 1e-12) / 10^5, nsmall = 5)
    }
    table_name <- paste0("qtable_c", name_search(c), "_lam", name_search(lambda))

    z <- get(table_name, envir = getNamespace("fdpbandsdata"))[
      d_max, paste0(1 - gamma)
    ]

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
      xi <- function(d) {
        if (d + 1 > d_max) {
          return(Inf)
        }
        max(0, floor(z * sqrt(B * (1 + B) * (d + 1)) + B * (d + 1) + 1e-9))
      }

      dbar <- function(threshold) {
        check_these <- c(
          1,
          setdiff(which(labels[1:threshold] == -1) - 1, 0),
          threshold
        )
        max_value <- 0
        for (t in check_these) {
          new_value <- ceiling(target_fun(t) - xi(decoy_fun(t)) - 1e-12)
          if (new_value > max_value) {
            max_value <- new_value
          }
        }
        return(max_value)
      }

      with_interpolation <- function(threshold) {
        if (threshold == 0) {
          0
        } else if (decoy_fun(threshold) >= 5e4) {
          1
        } else {
          (target_fun(threshold) - dbar(threshold)) /
            max(1, target_fun(threshold))
        }
      }

      return(vapply(thresholds, with_interpolation, numeric(1)))
    } else {
      without_interpolation <- function(d, t) {
        if (t == 0) {
          return(0)
        }
        if (d < 5e4) {
          min(1,
            floor(z * sqrt(B * (1 + B) * (d + 1)) + B * (d + 1) + 1e-12) /
              max(t, 1)
          )
        } else {
          1
        }
      }

      decoys <- vapply(thresholds, function(t) decoy_fun(t), numeric(1))
      targets <- vapply(thresholds, function(t) target_fun(t), numeric(1))

      return(mapply(without_interpolation, decoys, targets))
    }
  } else {
    stop(
      "The standardized band requires precomputed data tables. You may choose",
      " to run devtools::install_github(\"uni-Arya/fdpbandsdata\") to install",
      " these tables."
    )
  }
}

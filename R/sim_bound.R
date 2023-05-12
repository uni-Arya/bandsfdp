#' Simultaneous Band
#'
#' This function computes upper prediction bounds on the target wins among the
#' top \eqn{k} hypotheses of TDC, for each \eqn{k = 1,\ldots,n} where \eqn{n}
#' is the total number of hypotheses.
#'
#' @param labels A vector of (ordered) labels. See details below.
#' @param gamma The confidence parameter of the band. Typical values include
#' \code{gamma = 0.05} or \code{gamma = 0.01}.
#' @param type A character string specifying which band to use. Must be one of
#' \code{"stband"} or \code{"uniband"}.
#' @param d_max An optional positive integer specifying the maximum number
#' of decoy wins considered in calculating the bands.
#' @param max_fdp A number specifying the maximum FDP considered by the user in
#' calculating the bands. Used to compute \code{d_max} if \code{d_max} is
#' set to \code{NULL}.
#' @param c Determines the ranks of the target score that are considered
#' winning. Defaults to \code{c = 0.5} for (single-decoy) TDC.
#' @param lambda Determines the ranks of the target score that are
#' considered losing. Defaults to \code{lambda = 0.5} for (single-decoy) TDC.
#'
#' @details In (single-decoy) TDC, each hypothesis is associated to a
#' winning score and a label (1 for a target win, -1 for a decoy win). This
#' function assumes that the hypotheses are ordered in decreasing order of
#' winning scores (with ties broken at random). The argument \code{labels},
#' therefore, must be ordered according to this rule.
#'
#' This function also supports the extension of TDC that uses multiple
#' decoys. In that setup, the target score is competed with multiple decoy
#' scores and the rank of the target score after competition is used to determine whether the
#' hypothesis is a target win (label = 1), decoy win (-1) or uncounted (0).
#' The top \code{c} proportion of ranks are considered winning, the bottom
#' \code{1-lambda} losing, and all the rest uncounted.
#'
#' The threshold of TDC is given by the formula (assuming hypotheses are ordered):
#' \deqn{\max\{k : \frac{D_k + 1}{T_k \vee 1} \cdot \frac{c}{1-\lambda} \leq \alpha\}}{%
#' max{k : (D_k + 1)/max(T_k, 1) \le (1-\lambda)\alpha/c}}
#' where \eqn{T_k} is the number of target wins among the top
#' \eqn{k} hypotheses, and \eqn{D_k} is the number of decoy wins similarly.
#'
#' The argument \code{gamma} sets a confidence level of \code{1-gamma}. Both
#' the uniform and standardized bands require pre-computed Monte Carlo
#' statistics, so only certain values of \code{gamma} are available to use.
#' Commonly used confidence levels, like 0.95 and 0.99, are available.
#' We refer the reader to the README of this package for more details.
#'
#' The argument \code{d_max} controls the rate at which the returned bounds
#' increase: a larger \code{d_max} results in a more conservative bound.
#' If, however, \eqn{D_k + 1} exceeds \code{d_max} for some index \eqn{k}, each target
#' win thereafter is considered a false discovery when computing the bound.
#' Thus it is important that \code{d_max}, chosen a priori, is large enough. Given
#' it is sufficiently large, the precise value of \code{d_max} does not have a
#' significant effect on the resulting bounds (see <https://arxiv.org/abs/2302.11837> for more details).
#'
#' We recommend setting \code{d_max = NULL} so that it is computed automatically
#' using \code{max_fdp}. This argument ensures that \eqn{D_k + 1} never
#' exceeds \code{d_max} when the (non-interpolated) FDP bound on the top
#' \eqn{k} hypotheses is less than \code{max_fdp}.
#'
#' @return A vector of upper prediction bounds on the FDP of target wins among
#' the top \eqn{k} hypotheses for each \eqn{k = 1,\ldots,n} where \eqn{n}
#' is the total number of hypotheses.
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
#'   head(sim_bound(labels, gamma, "stband"))
#' }
#'
#' @references Ebadi et al. (2022), Bounding the FDP in competition-based
#' control of the FDR <https://arxiv.org/abs/2302.11837>.
sim_bound <- function(labels, gamma, type,
                      d_max = NULL, max_fdp = 0.5,
                      c = 0.5, lambda = 0.5) {
  if (!type %in% c("stband", "uniband")) {
    stop("Invalid type argument in sim_bound().")
  }
  if (requireNamespace("fdpbandsdata", quietly = TRUE)) {
    n <- length(labels)
    indices <- 1:n
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
        stop("Calculated d_max is 0. The bounds will all be 1.")
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

    decoys <- cumsum(labels == -1)
    targets <- cumsum(labels == 1)

    dbar <- function(type) {
      dbar_vec <- rep(0, n)

      # Compute d_bar for k = 1,...,n
      if (type == "stband") {
        for (k in 1:n) {
          if (k == 1) {
            dbar_vec[1] <- max(
              0,
              ceiling(targets[1] - xi_sb(decoys[1], d_max) - 1e-12)
            )
            next
          }
          dbar_vec[k] <- max(
            dbar_vec[k - 1],
            ceiling(targets[k] - xi_sb(decoys[k], d_max) - 1e-12)
          )
        }
      } else {
        for (k in 1:n) {
          if (k == 1) {
            dbar_vec[1] <- max(
              0,
              ceiling(targets[1] - xi_ub(decoys[1], d_max) - 1e-12)
            )
            next
          }
          dbar_vec[k] <- max(
            dbar_vec[k - 1],
            ceiling(targets[k] - xi_ub(decoys[k], d_max) - 1e-12)
          )
        }
      }

      return(dbar_vec)
    }

    with_interpolation <- function(indices, type) {
      bounds <- rep(0, length(indices))
      dbar_vec <- dbar(type)

      for (i in seq_along(indices)) {
        threshold <- indices[i]
        if (threshold == 0) {
          bounds[i] <- 0
        } else {
          bounds[i] <- (targets[i] - dbar_vec[i]) /
            max(1, targets[i])
        }
      }

      return(bounds)
    }

    return(with_interpolation(indices, type))
  } else {
    stop(
      "Simultaneous bands require precomputed data tables. You may choose",
      " to run devtools::install_github(\"uni-Arya/fdpbandsdata\") to install",
      " these tables."
    )
  }
}

#' @rdname sim_bound
simband <- sim_bound

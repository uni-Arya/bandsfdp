#' Uniform Band
#'
#' This function computes an upper prediction bound, derived from the
#' uniform band, on the FDP in TDC's list of discoveries.
#'
#' @param thresholds The rejection threshold of TDC. If given as
#' a vector, an upper prediction bound is returned for each element.
#' @param labels A vector of (ordered) labels. See details below.
#' @param alpha The FDR threshold.
#' @param gamma The confidence parameter of the bound. Typical values include
#' \code{gamma = 0.05} or \code{gamma = 0.01}.
#' @param c Determines the ranks of the target score that are considered
#' winning. Defaults to \code{c = 0.5} for (single-decoy) TDC.
#' @param lambda Determines the ranks of the target score that are
#' considered losing. Defaults to \code{lambda = 0.5} for (single-decoy) TDC.
#' @param n The number of hypotheses. Defaults to the length of \code{labels}.
#' @param interpolate A boolean indicating whether the bands should be
#' interpolated. Offers a slight boost in performance at the cost of computing
#' power. Defaults to \code{TRUE}.
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
#' The threshold of TDC is given by the formula:
#' \deqn{\max\{k : \frac{D_k + 1}{T_k \vee 1} \cdot \frac{c}{1-\lambda} \leq \alpha\}}{%
#' max{k : (D_k + 1)/max(T_k, 1) \le (1-\lambda)\alpha/c}}
#' where \eqn{T_k} is the number of target wins among the top
#' \eqn{k} hypotheses, and \eqn{D_k} is the number of decoy wins similarly.
#'
#' The argument \code{gamma} sets a confidence level of \code{1-gamma}. Since
#' the uniform band requires pre-computed Monte Carlo statistics, only
#' certain values of \code{gamma} are available to use. Commonly used
#' confidence levels, like 0.95 and 0.99, are available. We refer the reader
#' to the README of this package for more details.
#'
#' The argument \code{alpha}, used to compute the threshold of TDC, is also
#' used in this function. It serves to compute an appropriate d_max
#' for a non-trivial bound. In particular, if the user inputs a vector of
#' \code{thresholds}, a bound is returned for each element of
#' \code{thresholds} using the same d_max. For more details, see:
#' <https://arxiv.org/abs/2302.11837>.
#'
#' We recommend the use of \code{interpolate = TRUE} (default), as it generally
#' results in a tighter bound. This comes at the cost of performance: the bound
#' for each threshold is computed in O(n) time with interpolation and O(1)
#' without.
#'
#' @return An upper prediction bound on the FDP in TDC's list of discoveries.
#' If \code{thresholds} is a vector, returns an upper prediction bound for each
#' element of \code{thresholds}.
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
#'   uniband(thresholds, labels, alpha, gamma)
#' }
#'
#' @references Ebadi et al. (2022), Bounding the FDP in competition-based
#' control of the FDR <https://arxiv.org/abs/2302.11837>.
uniband <- function(thresholds, labels,
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

    # Look up the appropriate u_gamma table
    name_search <- function(num) {
      format(ceiling(num * 10^5 - 1e-12) / 10^5, nsmall = 5)
    }
    table_name <- paste0("utable_c", name_search(c), "_lam", name_search(lambda))

    rho <- get(table_name, envir = getNamespace("fdpbandsdata"))[[1]][
      d_max, paste0(1 - gamma)
    ]
    sigma <- get(table_name, envir = getNamespace("fdpbandsdata"))[[2]][
      d_max, paste0(1 - gamma)
    ]
    flip <- get(table_name, envir = getNamespace("fdpbandsdata"))[[3]][
      d_max, paste0(1 - gamma)
    ]
    u <- ifelse(stats::runif(1) < flip, rho, sigma)

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
        floor(stats::qnbinom(1 - u, d + 1, prob = (1 / (1 + B))) + 1e-9)
      }

      dbar <- function(threshold) {
        check_these <- c(
          1,
          setdiff(which(labels[1:threshold] == -1) - 1, 0),
          threshold
        )
        max_value <- 0
        for (t in check_these) {
          new_value <- ceiling(target_fun(t) - xi(decoy_fun(t))  - 1e-12)
          if (new_value > max_value) {
            max_value <- new_value
          }
        }
        return(max_value)
      }

      with_interpolation <- function(threshold) {
        if (threshold == 0) {
          0
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
            floor(stats::qnbinom(1 - u, d + 1, prob = (1 / (1 + B))) + 1e-12) /
              max(1, t)
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
      "The uniform band requires precomputed data tables. You may choose",
      " to run devtools::install_github(\"uni-Arya/fdpbandsdata\") to install",
      " these tables."
    )
  }
}

#' @rdname uniband
tdc_ub <- uniband

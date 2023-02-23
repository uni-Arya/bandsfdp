#' Calculate the band of Katsevich-Ramdas on TDC's list of discoveries
#'
#' In TDC, each hypothesis is attached to a winning score and a label (1 for a
#' target win and -1 for a decoy win). \code{krband()} assumes that the
#' hypotheses are ordered in decreasing order of winning scores (breaking
#' ties at random). In particular, the argument \code{labels} is ordered
#' according to this rule.
#'
#' @param thresholds A vector of rejection thresholds resulting from TDC.
#' @param labels A vector of (ordered) labels resulting from TDC.
#' @param gamma The confidence level of the band. Typical values include
#' \code{gamma = 0.05} or \code{gamma = 0.01}.
#' @param c Determines the ranks of the target score that are considered
#' winning. Defaults to \code{c = 0.5} for TDC.
#' @param lambda Determines the ranks of the target score that are
#' considered losing. Defaults to \code{lambda = 0.5} for TDC.
#' @param interpolate A boolean indicating whether the bands should be
#' interpolated. Offers a slight boost in performance at the cost of computing
#' power. Defaults to \code{TRUE}.
#'
#' @return A vector of upper prediction bounds on the FDP for each rejection
#' threshold.
#' @export
#'
#' @examples
#' set.seed(123)
#' thresholds <- c(250, 500, 750, 1000)
#' labels <- c(
#'   rep(1, 250),
#'   sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.9, 0.1)),
#'   sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.5, 0.5)),
#'   sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.1, 0.9))
#' )
#' gamma <- 0.05
#' krband(thresholds, labels, gamma)
krband <- function(thresholds, labels,
                   gamma,
                   c = 0.5, lambda = 0.5,
                   interpolate = TRUE) {
  if (any(thresholds < 0)) {
    stop("Negative threshold detected.")
  }
  B <- c / (1 - lambda)
  C <- -log(gamma) / log(1 + (1 - gamma^B) / B)

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
      floor(C * (1 + B * d) + 1e-9)
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
        return(0)
      }
      (target_fun(threshold) - dbar(threshold)) / max(1, target_fun(threshold))
    }

    return(vapply(thresholds, with_interpolation, numeric(1)))
  } else {
    without_interpolation <- function(d, t) {
      if (t == 0) {
        return(0)
      }
      min(1, floor(C * (1 + B * d) + 1e-12) / max(1, t))
    }

    decoys <- vapply(thresholds, function(t) decoy_fun(t), numeric(1))
    targets <- vapply(thresholds, function(t) target_fun(t), numeric(1))

    return(mapply(without_interpolation, decoys, targets))
  }
}

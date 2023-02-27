
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bandsfdp

<!-- badges: start -->
<!-- badges: end -->

This package provides several different upper prediction bands on the
FDP in the list of discoveries resulting from competition-based setups
(see [Ebadi et al. (2022)](https://arxiv.org/abs/2302.11837)). Such
setups include target-decoy competition (TDC) in computational mass
spectrometry and the knockoff construction in regression.

The band provides an upper prediction bound on the FDP: given a
rejection threshold returned by TDC, a list of target/decoy labels, and
a desired confidence level
![1 - \gamma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%20-%20%5Cgamma "1 - \gamma"),
these functions return a real number
![\eta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ceta "\eta")
such that the FDP in the list of discoveries is
![\leq \eta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleq%20%5Ceta "\leq \eta")
with probability
![\geq 1 - \gamma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cgeq%201%20-%20%5Cgamma "\geq 1 - \gamma").

## Installation

You can install the development version of bandsfdp from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("uni-Arya/bandsfdp")
```

## Usage

The following inputs, pertaining to the competition-based procedure, are
required:

1.  One, or several, rejection `thresholds` (the user may choose their
    own threshold(s) as well).
2.  A list of `labels`
    (![-1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-1 "-1")
    for a decoy win,
    ![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1 "1")
    for a target win) that are ordered so the corresponding winning
    scores are decreasing.

Note to compute the standardized band (`stband()`) and uniform band
(`uniband()`), the user must download a set of pre-computed Monte-Carlo
data tables using `devtools::install_github("uni-Arya/fdpbandsdata")`.

### With a Single Decoy Score (as in TDC)

By default, each band assumes the use of a single decoy score, as in
TDC. In this case, in addition to the `thresholds` and `labels`, one
must specify the confidence parameter `gamma` (for a `1 - gamma`
confidence level). Note that the functions currently support
`gamma = 0.01, 0.025, 0.5, 0.1, 0.8, 0.5`, but more data can be
generated using the source code available at
[GitHub](https://github.com/uni-Arya/fdpbandsdata). For the standardized
and uniform band, the value of `alpha` (the FDR tolerance) must also be
specified.

``` r
suppressPackageStartupMessages(library(bandsfdp))
set.seed(123)

if (requireNamespace("fdpbandsdata", quietly = TRUE)) {
  thresholds <- c(250, 500, 750, 1000)
  labels <- c(
    rep(1, 250),
    sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.9, 0.1)),
    sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.5, 0.5)),
    sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.1, 0.9))
  )
  alpha <- 0.05
  gamma <- 0.05
  
  print(stband(thresholds, labels, alpha, gamma))
  print(uniband(thresholds, labels, alpha, gamma))
  print(krband(thresholds, labels, gamma))
}
#> [1] 0.02000000 0.09453782 0.26825127 0.29575163
#> [1] 0.02400000 0.08823529 0.26315789 0.29084967
#> [1] 0.0160000 0.2184874 0.3684211 0.3921569
```

### With Multiple Decoy Scores

For competition-based setups that utilizes multiple decoys, one must
specify two additional parameters
![c \leq \lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c%20%5Cleq%20%5Clambda "c \leq \lambda")
of the form
![k/(d + 1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%2F%28d%20%2B%201%29 "k/(d + 1)")
where
![d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d "d")
is the number of decoy scores used for each hypothesis and
![1 \leq k \leq d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%20%5Cleq%20k%20%5Cleq%20d "1 \leq k \leq d")
is an integer.

For example, if we use
![3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;3 "3")
decoy scores for each hypothesis, we may take
![c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c "c")
and
![\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
to be
![1/4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%2F4 "1/4"),
![1/2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%2F2 "1/2"),
or
![3/4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;3%2F4 "3/4"),
subject to
![c \leq \lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c%20%5Cleq%20%5Clambda "c \leq \lambda").
Briefly, the value of
![c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c "c")
determines the ranks in which the target score is considered “winning.”
E.g., if
![c = 1/4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c%20%3D%201%2F4 "c = 1/4"),
![H_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;H_i "H_i")
is labelled as a target win whenever its corresponding target score is
the highest ranked score among all targets/decoys for that hypothesis.
Similarly,
![\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
determines when the target score is “losing.” For more details, see
[Ebadi et al. (2022)](https://arxiv.org/abs/2302.11837).

``` r
suppressPackageStartupMessages(library(bandsfdp))
set.seed(123)

if (requireNamespace("fdpbandsdata", quietly = TRUE)) {
  thresholds <- c(250, 500, 750, 1000)
  labels <- c(
    rep(1, 250),
    sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.9, 0.1)),
    sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.5, 0.5)),
    sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.1, 0.9))
  )
  alpha <- 0.05
  gamma <- 0.05
  c <- 0.25
  lambda <- 0.25
  
  print(stband(thresholds, labels, alpha, gamma, c, lambda))
  print(uniband(thresholds, labels, alpha, gamma, c, lambda))
  print(krband(thresholds, labels, gamma, c, lambda))
}
#> [1] 0.00800000 0.03991597 0.16298812 0.19444444
#> [1] 0.01200000 0.03781513 0.15449915 0.18627451
#> [1] 0.00800000 0.05042017 0.22750424 0.25653595
```

### Interpolated Bands

By default, all bands are interpolated, which requires the computation
of a running maximum. This generally results in a slightly tighter band,
but at the cost of computational power. We recommend the use of
`interpolate = TRUE`, unless it is too time-consuming. The code below
shows an example using non-interpolated bands.

``` r
suppressPackageStartupMessages(library(bandsfdp))
set.seed(123)

if (requireNamespace("fdpbandsdata", quietly = TRUE)) {
  thresholds <- c(250, 500, 750, 1000)
  labels <- c(
    rep(1, 250),
    sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.9, 0.1)),
    sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.5, 0.5)),
    sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.1, 0.9))
  )
  alpha <- 0.05
  gamma <- 0.05
  c <- 0.25
  lambda <- 0.25
  
  print(stband(thresholds, labels, alpha, gamma, c, lambda, interpolate = FALSE))
  print(uniband(thresholds, labels, alpha, gamma, c, lambda, interpolate = FALSE))
  print(krband(thresholds, labels, gamma, c, lambda, interpolate = FALSE))
}
#> [1] 0.00800000 0.03991597 0.13921902 0.28267974
#> [1] 0.01200000 0.03781513 0.12903226 0.26633987
#> [1] 0.00800000 0.05252101 0.26146010 0.59967320
```

### Simultaneous FDP bounds

Rather than compute a bound on the FDP in TDC’s list of discoveries, one
may be interested in computing a bound on the FDP on the top
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")
hypotheses for all
![i = 1, \ldots, n =: \[n\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%20%3D%201%2C%20%5Cldots%2C%20n%20%3D%3A%20%5Bn%5D "i = 1, \ldots, n =: [n]"),
where
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
is the number of hypotheses (or for
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")
in a subset of
![\[n\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Bn%5D "[n]")).
In this case, the function `simband()` may be used. This requires the
following arguments:

- The set of (ordered) `labels`, confidence parameter `gamma`,
  parameters `c` and `lambda`, and a boolean `interpolate`, as described
  in the previous sections.
- A character argument `type` which is either `"stband"` or `"uniband"`,
  specifying the type of band to be used to compute simultaneous FDP
  bounds.
- The set of desired `indices` (denoted by
  ![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")
  above).
- The maximum number of decoy wins considered for the bands `d_max`
  (defaults to `NULL`, in which case is automatically computed using
  `max_fdp` below).
- The maximum considered FDP for the simultaneous bands `max_fdp`
  (defaults to `max_fdp = 0.5`).

The arguments `d_max` and `max_fdp` control the rate at which the
simultaneous bounds are increasing. We refer the reader to Section 3 of
[Ebadi et al. (2022)](https://arxiv.org/abs/2302.11837) for more
details. Note that `max_fdp` is the
![\alpha](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha "\alpha")
used to calculate
![d\_\infty](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d_%5Cinfty "d_\infty").

``` r
suppressPackageStartupMessages(library(bandsfdp))
set.seed(123)

if (requireNamespace("fdpbandsdata", quietly = TRUE)) {
  set.seed(123)
  labels <- c(
    rep(1, 250),
    sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.9, 0.1)),
    sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.5, 0.5)),
    sample(c(1, -1), size = 250, replace = TRUE, prob = c(0.1, 0.9))
  )
  gamma <- 0.05
  print(simband(labels, gamma, "stband", indices = c(100, 250, 500, 1000)))
}
#> [1] 0.05000000 0.02000000 0.09663866 0.29738562
```

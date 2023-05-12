# bandsfdp 1.1.0

## New features

- Included function `gen_bound`, which computes an upper prediction bound on the FDP on any arbitrary set of hypotheses $R$.
- Updated README to include an example of use.

## Minor improvements and fixes

- Updated aliases of functions so that `tdc_ub`/`tdc_sb`/`sim_bound` are the primary aliases, instead of `uniband`/`stband`/`simband`. Note both names can still be used to call these functions.
- Fixed a bug with `tdc_ub` and `tdc_sb` that returned a non-trivial bound when the number of decoys was larger than `d_max`

# bandsfdp 1.0.0

-   Initial release of the bandsfdp package.

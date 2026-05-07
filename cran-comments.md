## Resubmission spBPS (version 2.0-1):

-   Major release (previous was version 1.0-1).

-   New features:

    -   `spBPS()` now replaces the previous implementation with a fully optimised version based on OpenMP parallelism and native C++ solvers (no CVXR dependency for stacking weight estimation).

    -   `predict.spBPS()`: new S3 method for posterior predictive sampling at new spatial locations. Uses a streaming algorithm that never allocates u x u matrices, supporting arbitrarily large prediction sets with controlled memory usage via `pred_batch_size`.

-   Bug fixes and correctness:

    -   Fixed parameterisation of the inverse-Wishart sampler.

    -   Fixed matrix-t log-density evaluation to match the Sylvester identity formulation (no mniw dependency for matrix distribution sampling).

### R CMD check results

0 errors \| 0 warnings \| 1 note

-   NOTE: (checking for future file timestamps ... NOTE / unable to verify current time)

## Resubmission spBPS (version 1.0-1):

-   Major release (previous was version 0.0-4)

-   Addressed all the CRAN reviewer (Kurt Hornik) recommendations:

    -   Resolved R-devel check Issues.

-   Fixed compatibility with CVXR \>= 1.8 (psolve migration).

-   Introduced new functions:

    -   spBPS(): unified orchestrator function to simplify workflow execution.

-   Updated examples and vignette.

### R CMD check results

0 errors \| 0 warnings \| 1 note

-   NOTE: (checking for future file timestamps ... NOTE / unable to verify current time)

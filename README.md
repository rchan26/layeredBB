# layeredBB

Code to implement methods for simulating Brownian motion, Brownian bridges and layered Brownian bridge skeletons.
Please see function documentation for example usage.

## Installation

Simply run: `devtools::install_github('rchan26/layeredBB')`.

Note that there are vignettes available, but by default `devtools::install_github()` does not build them since they are time consuming. To force building, use `devtools::install_github('rchan26/layeredBB', build_vignettes = TRUE)`.

## Using with `parallel`

Here is an example to use in a cluster using [`parallel`](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf):

```
library(layeredBB)
library(parallel)

# creating parallel cluster with maximum number of cores available
n_cores <- parallel::detectCores()
cl <- parallel::makeCluster(n_cores)

# export functions from layeredBB into the cluster environment
parallel::clusterExport(cl, varlist = ls("package:layeredBB"))

##### do some stuff in cluster using layeredBB functions ######

# stop cluster
parallel::stopCluster(cl)
```

## Development workflow

If any code has been modified, call:

```
Rcpp::compileAttributes()
pkgbuild::compile_dll()
devtools::document()
devtools::install()
```

To build vignettes use:
```
devtools::build_vignettes()
```

Or, you can use
```
devtools::build()
```
to create a package bundle with the vignettes included.

Then can test the package by calling: `devtools::test()`

## Resources

* [On the exact and epsilon-strong simulation of (jump) diffusions](https://warwick.ac.uk/fac/sci/statistics/staff/academic-research/johansen/publications/PJR16.pdf)
* [Some Monte Carlo methods for jump diffusions](http://wrap.warwick.ac.uk/60602/)

## License

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

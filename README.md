## To install, use `install_github` from `devtools` package

```
devtools::install_github('rchan26/layeredBB')
```

## Using with `parallel`

Here is an example to use in a cluster using `parallel`:

```
  # load in package
  library(layeredBB)

  # creating parallel cluster with maximum number of cores
  n_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(n_cores)
  
  # export functions from layeredBB into the cluster environment
  parallel::clusterExport(cl, varlist = ls("package:layeredBB"))
  
  ##### do some stuff in cluster using layeredBB functions ######
  
  # stop cluster
  parallel::stopCluster(cl)
```


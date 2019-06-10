# To load into R: run `layers_rcpp.R`

```
devtools::install_github('rchan26/layeredBB')
```

Must have Rcpp package and a C++ compiler. This package uses the C compiler ("cc" or "gcc") and was designed using for Unix/Linux/Mac machines, where it should work without difficulty. With Mac OS X, it is necessary to first install `gcc`. 

I'm not sure if this package runs on Microsoft Windows. 

## Warning: using with `parallel`

If you want to use this with the `parallel` package, you must download the package and install it. To use it then must call `library(layeredBB)` - see [here](https://stackoverflow.com/questions/6074310/using-rcpp-within-parallel-code-via-snow-to-make-a-cluster).

Here is an example to use in a cluser:

```
  # load in package
  library(layeredBB)

  # creating parallel cluster
  n_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(n_cores)
  
  # export functions from layeredBB into the cluster environment
  parallel::clusterExport(cl, varlist = ls("package:layeredBB"))
  
  ##### do some stuff in cluster ######
  
  # stop cluster
  parallel::stopCluster(cl)
```

## Steps to install the package locally

#### On the command line (move to directory containing the package):
```
R CMD INSTALL layeredBB_1.0.tar.gz 
```

#### In R/RStudio:
```
library(devtools)
library(Rcpp)
Rcpp::compileAttributes('layeredBB')
devtools::build('layeredBB')
devtools::install('layeredBB')
library(layeredBB)
```

## Steps to create a package again on the command line

### Step One: create a placeholder package for Rcpp
```
Rscript -e "Rcpp::Rcpp.package.skeleton('layeredBB/')"
```

### Step Two: copy your c++ source and header files into the package folder

### Step Three: run 'compileAttributes()' to update the exports
```
Rscript -e "Rcpp::compileAttributes('layeredBB/')"
```

### Step Four: Build the package
```
R CMD build layeredBB/
```

### Step Five: Install the package
```
R CMD INSTALL layeredBB_1.0.tar.gz 
```

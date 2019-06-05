# To load into R: run layers_rcpp.R

```
source('https://raw.githubusercontent.com/rchan26/bessel_layers_rcpp/master/layers_rcpp.R')
```

This downloads layers_rcpp.cpp from this page and runs sourceCpp.

Must have Rcpp package and a C++ compiler. This package uses the C compiler ("cc" or "gcc") and was designed using for Unix/Linux/Mac machines, where it should work without difficulty. With Mac OS X, it is necessary to first install gcc. 

I'm not sure if this package runs on Microsoft Windows.

## Steps to install the package

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

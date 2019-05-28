# Steps to install the package

#### On the command line (move to directory containing the package):
```
R CMD INSTALL layeredBB_1.0.tar.gz 
```

#### In R/RStudio:
```
library(devtools)
devtools::install('layeredBB_1.0.tar.gz')
```

# Steps to create a package again on the command line

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



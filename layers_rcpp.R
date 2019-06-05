download_cpp <- function(url) { 
  tf <- tempfile(fileext = '.cpp') 
  download.file(url, tf, quiet=TRUE) 
  return(tf)
}

library(Rcpp)
url <- "https://raw.githubusercontent.com/rchan26/bessel_layers_rcpp/master/layers_rcpp.cpp"
sourceCpp(download_cpp(url))

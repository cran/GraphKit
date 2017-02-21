#' GraphKit
#' 
#'Efficient methods for constructing confidence intervals of monotone 
#'graph invariants, as well as testing for monotone graph properties. 
#'Many packages are available to estimate precision matrices, this 
#'package serves as a tool to extract structural properties from 
#'their induced graphs. By iteratively bootstrapping on only the 
#'relevant edge set, we are able to obtain the optimal interval size.
#' 
#' @docType package
#' @author Manbir Gulati <manbirgulati@gmail.com>, Junwei Lu, Han Liu
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib GraphKit, .registration=TRUE
#' @name GraphKit
#'
NULL


library(arm)
library(foreign)
library(bbmle)
library(testthat)
library(tictoc)
library(rstan)
library(jsonlite)
library(parallel)
library(rethinking) # devtools::install_github("rmcelreath/rethinking")
library(rake) # devtools::install_github("babeheim/rake")

n_cores <- 3

n_iter <- 4000
n_warmup <- floor(n_iter / 2)

beta_tol <- 0.05

init <- "random"

randomize_ids <- function(ids, reserved, seed = NA, nchars = 5) {
  if(all(is.na(ids)) | length(ids) == 0) stop("not valid ids")
  ids <- as.character(ids)
  id_list <- sort(unique(ids))
  new_ids <- id_maker(n = length(id_list), seed = seed, nchars = nchars)
  out <- new_ids[match(ids, id_list)]
  expect_identical(table(table(ids)), table(table(out)))
  expect_identical(which(is.na(ids)), which(is.na(out)))
  return(out)
}

id_maker <- function(n, reserved = "", seed = NA, nchars = NA){
  my_let <- letters 
  my_num <- 0:9 
  if(is.na(seed) | !is.numeric(seed)) set.seed(as.numeric(as.POSIXlt(Sys.time())))
  if(!is.na(seed) & is.numeric(seed)) set.seed(seed)
  output <- replicate(n, paste(sample(c(my_let, my_num), nchars, replace=TRUE), 
    collapse=''))
  rejected <- duplicated(output) | output %in% reserved |
    substr(output, 1, 1) %in% my_num
  while (any(rejected)) {
    output <- output[-which(rejected)]
    remaining <- n - length(output)
    output <- c(output, replicate(remaining, paste(sample(c(my_let, my_num), nchars, 
      replace=TRUE), collapse="")))
    rejected <- duplicated(output) | output %in% reserved |
      substr(output, 1, 1) %in% my_num
  }
  output
}

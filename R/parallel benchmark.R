library(doParallel)
library(doMC)
library(rbenchmark)
library(ggplot2)
cl <- makeCluster(4)
registerDoParallel(cl)
getDoParWorkers()
set.seed(1990)

####################
########## simple functions 

FUN <- function(x) { round(sqrt(x), 4) }
a <- lapply(1:10, function(i) i)

test1 <- benchmark("lapply" = lapply(1:10, FUN = FUN), 
                   "For loop" = for(i in 1:10){ FUN(i)},
                   "Foreach dopar" = foreach(i = 1:10) %dopar% FUN(i),
                   "Foreach do" = foreach(i = 1:10) %do% FUN(i),
                   "parLapply" = parLapply(cl = cl, X = a, fun = FUN),
                   "parSapply" = parSapply(cl = cl, X = a, FUN = FUN),
                   columns=c('test', 'elapsed', 'replications'),
                   replications = c(100, 200, 500, 1000))

ggplot() +
  geom_line(aes(x = replications, y = elapsed, colour = test), data = test1)

######################
########## matrix multiplication

### element wise multiplication
FUN <- function(x){
  a <- matrix(rnorm(10000), 100, 100)
  b <- matrix(rnorm(10000), 100, 100)
  a * b
}

test2 <- benchmark("lapply" = lapply(1:10, FUN = FUN), 
                   "For loop" = for(i in 1:10){ FUN(i)},
                   "Foreach dopar" = foreach(i = 1:10) %dopar% FUN(i),
                   "Foreach do" = foreach(i = 1:10) %do% FUN(i),
                   "parLapply" = parLapply(cl = cl, X = a, fun = FUN),
                   "parSapply" = parSapply(cl = cl, X = a, FUN = FUN),
                   columns=c('test', 'elapsed', 'replications'),
                   replications = c(100, 200, 500, 1000))

ggplot() +
  geom_line(aes(x = replications, y = elapsed, colour = test), data = test2)



### inner product
FUN <- function(x) {
  a <- matrix(rnorm(10000), 100, 100)
  b <- matrix(rnorm(10000), 100, 100)
  a %*% b
}

test3 <- benchmark("lapply" = lapply(1:10, FUN = FUN), 
                   "For loop" = for(i in 1:10){ FUN(i)},
                   "Foreach dopar" = foreach(i = 1:10) %dopar% FUN(i),
                   "Foreach do" = foreach(i = 1:10) %do% FUN(i),
                   "parLapply" = parLapply(cl = cl, X = a, fun = FUN),
                   "parSapply" = parSapply(cl = cl, X = a, FUN = FUN),
                   columns=c('test', 'elapsed', 'replications'),
                   replications = c(100, 200, 500, 1000))

ggplot() +
  geom_line(aes(x = replications, y = elapsed, colour = test), data = test3)


#### inner product with unbalanced dimensions
FUN <- function(x) {
  a <- matrix(rnorm(10000), 1000, 10)
  b <- matrix(rnorm(10000), 10, 1000)
  a %*% b
}

test4 <- benchmark("lapply" = lapply(1:10, FUN = FUN), 
                   "For loop" = for(i in 1:10){ FUN(i)},
                   "Foreach dopar" = foreach(i = 1:10) %dopar% FUN(i),
                   "Foreach do" = foreach(i = 1:10) %do% FUN(i),
                   "parLapply" = parLapply(cl = cl, X = a, fun = FUN),
                   "parSapply" = parSapply(cl = cl, X = a, FUN = FUN),
                   columns=c('test', 'elapsed', 'replications'),
                   replications = c(10, 20, 50, 100))

ggplot() +
  geom_line(aes(x = replications, y = elapsed, colour = test), data = test4)

###### glm
FUN <- function(i) {
  ind <- sample(100, 100, replace=TRUE)
  result1 <- glm(Species~Sepal.Length, family=binomial(logit), data = iris[ind,])
  coefficients(result1)
}

test5 <- benchmark("lapply" = lapply(1:10, FUN = FUN), 
                   "For loop" = for(i in 1:10){ FUN(i)},
                   "Foreach dopar" = foreach(i = 1:10) %dopar% FUN(i),
                   "Foreach do" = foreach(i = 1:10) %do% FUN(i),
                   "parLapply" = parLapply(cl = cl, X = a, fun = FUN),
                   "parSapply" = parSapply(cl = cl, X = a, FUN = FUN),
                   columns=c('test', 'elapsed', 'replications'),
                   replications = c(100, 200, 500, 100))

ggplot() +
  geom_line(aes(x = replications, y = elapsed, colour = test), data = test5)

###### glm but in a 'list' fashion
a <- lapply(1:10, function(i) iris)

FUN <- function(x) {
  ind <- sample(100, 100, replace=TRUE)
  result1 <- glm(Species~Sepal.Length, family=binomial(logit), data = x[ind,])
  coefficients(result1)
}

test6 <- benchmark("lapply" = lapply(a, FUN = FUN), 
                   "For loop" = for(i in 1:10){ FUN(a[[i]])},
                   "Foreach dopar" = foreach(i = 1:10) %dopar% FUN(a[[i]]),
                   "Foreach do" = foreach(i = 1:10) %do% FUN(a[[i]]),
                   "parLapply" = parLapply(cl = cl, X = a, fun = FUN),
                   "parSapply" = parSapply(cl = cl, X = a, FUN = FUN),
                   columns=c('test', 'elapsed', 'replications'),
                   replications = c(100, 200, 500, 1000))

ggplot() +
  geom_line(aes(x = replications, y = elapsed, colour = test), data = test6)

####### compare by using doMC as backend
registerDoMC(cores = 4)

test7 <- benchmark("lapply" = lapply(a, FUN = FUN), 
                   "For loop" = for(i in 1:10){ FUN(a[[i]])},
                   "Foreach dopar" = foreach(i = 1:10) %dopar% FUN(a[[i]]),
                   "Foreach do" = foreach(i = 1:10) %do% FUN(a[[i]]),
                   "parLapply" = parLapply(cl = cl, X = a, fun = FUN),
                   "parSapply" = parSapply(cl = cl, X = a, FUN = FUN),
                   columns=c('test', 'elapsed', 'replications'),
                   replications = c(100, 200, 500, 1000))

ggplot() +
  geom_line(aes(x = replications, y = elapsed, colour = test), data = test7)



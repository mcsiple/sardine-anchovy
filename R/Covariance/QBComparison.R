# Compare influence of Q and B matrices on covariance

vec <- function(x) { #converts matrix to vector
  if(class(x) != "matrix") {
    stop("`x` must be a matrix")
  }
  dim(x) <- c(prod(dim(x)), 1)
  return(x)
}

ivec <- function(x, cols = NULL) { # turns vector into matrix
  if(is.null(cols)) {
    rows <- cols <- sqrt(length(x))
    if(rows != round(rows)) {
      stop("The argument `cols` is undefined and the length of `x` does not have an integer square root.")
    }
  } else {
    rows <- length(x) / cols
    if(rows != round(rows)) {
      stop("The quantity `rows = length(x) / cols` is not an integer.")
    }
  }
  dim(x) <- c(rows, cols)
  return(x)
}


BB <- matrix(c(0.98, 0.03, -0.1, 0.67), 2, 2) #global B
BB <- diag(bb)


QQ_Benguela <- matrix(c(0.08, 0.06, 0.06, 0.37), 2, 2) # process var/cov
# [,1] [,2]
# [1,] 0.08 0.06
# [2,] 0.06 0.37
QQ_California <- matrix(c(0.12, -0.11, -0.11, 0.48), 2, 2) #estimates from step 1, MARSS
QQ_Humboldt <- matrix(c(0.19, -0.03, -0.03, 0.37), 2, 2)
QQ_NEAtlantic <- matrix(c(0.11, -0.15, -0.15, 0.49), 2, 2)
QQ_KuroshioOyashio <- matrix(c(0.11, 0.10, 0.10, 0.49), 2, 2)

# Set these differently depending on which system you're interested in
#QQ <- QQ_California
#QQ <- QQ_Humboldt
#QQ <- QQ_NEAtlantic
QQ <- QQ_KuroshioOyashio

II <- diag(2)
##      [,1] [,2]
## [1,]    1    0
## [2,]    0    1

vec_SS <- solve(II %x% II - BB %x% BB) %*% vec(QQ)
##            [,1]
## [1,]  0.9803922
## [2,] -0.8722110
## [3,] -0.8722110
## [4,]  1.4097363

SS <- ivec(vec_SS)
##            [,1]      [,2]
## [1,]  0.9803922 -0.872211
## [2,] -0.8722110  1.409736

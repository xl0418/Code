noplotltt95 <- function (trees, alpha = 0.05, log = FALSE, method = c("lineages", 
                                                       "times"), mode = c("median", "mean"), ...) 
{
  if (!inherits(trees, "multiPhylo")) 
    stop("trees should be an object of class \"multiPhylo\".")
  method <- method[1]
  mode <- mode[1]
  if (hasArg(res)) 
    res <- list(...)$res
  else res <- 100
  X <- ltt(trees, plot = FALSE, gamma = FALSE)
  if (method == "times") {
    N <- length(X)
    tt <- sapply(X, function(x) max(x$times))
    zz <- max(tt) - tt
    for (i in 1:N) X[[i]]$times <- X[[i]]$times + zz[i]
    n <- sapply(X, function(x) max(x$ltt))
    if (all(n == max(n))) 
      n <- max(n)
    else stop("for method=\"times\" all trees must contain the same numer of lineages")
    LL <- sapply(X, function(x) x$times[1:length(x$times)])
    ii <- floor(alpha/2 * N)
    jj <- ceiling((1 - alpha/2) * N)
    low <- apply(LL, 1, function(x) sort(x)[ii])
    high <- apply(LL, 1, function(x) sort(x)[jj])
    ll <- if (mode == "median") 
      apply(LL, 1, function(x) median(x)[1])
    else rowMeans(LL)
    obj <- cbind(c(1:n, n), low, ll, high)
    colnames(obj) <- c("lineages", "low(time)", "time", 
                       "high(time)")
    rownames(obj) <- NULL
  }
  else if (method == "lineages") {
    N <- length(X)
    tt <- sapply(X, function(x) max(x$times))
    zz <- max(tt) - tt
    for (i in 1:N) X[[i]]$times <- X[[i]]$times + zz[i]
    tt <- 0:res * max(tt)/res
    ll <- low <- high <- vector()
    for (i in 1:length(tt)) {
      ss <- vector()
      for (j in 1:N) {
        ii <- 2
        while (tt[i] > X[[j]]$times[ii] && ii < length(X[[j]]$times)) ii <- ii + 
            1
        ss[j] <- X[[j]]$ltt[ii - 1]
      }
      ll[i] <- if (mode == "median") 
        median(ss)
      else mean(ss)
      low[i] <- sort(ss)[floor(alpha/2 * N)]
      high[i] <- sort(ss)[ceiling((1 - alpha/2) * N)]
    }
    obj <- cbind(tt, low, ll, high)
    colnames(obj) <- c("time", "low(lineages)", "lineages", 
                       "high(lineages)")
    rownames(obj) <- NULL
  }
  attr(obj, "class") <- "ltt95"
  attr(obj, "alpha") <- alpha
  attr(obj, "method") <- method
  attr(obj, "mode") <- mode
  attr(obj, "log") <- log
  invisible(obj)
}

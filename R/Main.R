#' Correlation between rows in x vs. y
#' @param x matrix whose rows are to be correlated with y
#' @param y vector of numbers
#' @param ... extra arguments to pass to cor function
#' @return correlation between each row and y
#' @importFrom stats cor
#' @export
cor2 <- function(x, y, ...) {
  temp <- as.numeric(stats::cor(t(x), y, ...))
  names(temp) <- rownames(x)
  temp
}

#' Correlation between rows in x vs. y
#' @param x matrix whose rows to be correlated with y
#' @param y vector of numbers
#' @param n number of trials for confidence intervals and P values
#' @param ... extra arguments to pass to cor function
#' @return data frame of correlations, confidence intervals, P value, false discovery rate
#' @examples
#' cor_test(t(iris[,1:3]), iris[,4])
#' @importFrom boot boot
#' @importFrom stats p.adjust
#' @export
cor_test <- function(x, y, n=1000, ...) {
  if (is.vector(x)) {
    # Change into matrix
    x <- t(matrix(x))
  }
  est <- cor2(x, y, ...)
  f <- function(d, ind) {
    cor2(t(d[ind,]), y[ind], ...)
  }
  b.obj <- boot::boot(t(x), f, R=n)
  # Get confidence intervals from bootstraps
  ci <- do.call(rbind, lapply(1:nrow(x), function(i) {
    boot::boot.ci(b.obj, type="perc", index=i)$percent[,4:5]
  }))
  # Permutation test to get P value
  p <- perm_test(function(y2) cor2(x, y2, ...), y, function() sample(y), n=n)
  # Show the results
  data.frame(
    Estimate = est,
    Lower = ci[,1],
    Upper = ci[,2],
    P = p,
    FDR = stats::p.adjust(p, "fdr")
  )
}

#' Fits y to x using robust linear regression and returns regression coefficients
#' @param x matrix of values where rows = genes, columns = samples
#' @param y value for every sample in x
#' @return coefficients from robust linear regression
#' @importFrom MASS rlm
#' @importFrom stats coef
#' @examples
#' rlm_coef(t(iris[,1:3]), iris[,4])
#' @export
rlm_coef <- function(x, y) {
  if (ncol(x) != length(y)) {
    stop("There must be as many y values as number of samples.")
  }
  model <- MASS::rlm(matrix(y) ~ t(x))
  coef.vec <- stats::coef(model)
  names(coef.vec) <- gsub("^t\\(x\\)", "", names(coef.vec))
  coef.vec
}

#' Fits robust linear model on all variables and
#'     tests whether coefficients of terms in model are equal to 0
#' @param x matrix of numbers whose rows are to be regressed with y
#' @param y vector of numbers
#' @param n number of trials for confidence interval and P values
#' @return coefficients, confidence intervals, and P values from linear model fit
#' @examples
#' rlm_coef_test(t(iris[,1:3]), iris[,4])
#' @importFrom boot boot boot.ci
#' @importFrom stats p.adjust
#' @export
rlm_coef_test <- function(x, y, n=1000) {
  # Bootstrap confidence intervals
  f <- function(d, ind) {
    rlm_coef(t(d[ind,]), y[ind])
  }
  b.obj <- boot::boot(t(x), f, R=n)
  # Get confidence intervals from bootstraps
  ci <- do.call(rbind, lapply(1:(nrow(x)+1), function(i) {
    boot::boot.ci(b.obj, type="perc", index=i)$percent[,4:5]
  }))
  # Prepare output
  est <- rlm_coef(x, y)
  lower <- ci[,1]
  upper <- ci[,2]
  # Permutation test by scrambling y
  p <- perm_test(function(y2) rlm_coef(x, y2), y, function() sample(y), n=n)
  fdr <- stats::p.adjust(p, "fdr")
  temp <- data.frame(
    Estimate=est,
    Lower=lower,
    Upper=upper,
    P=p,
    FDR=fdr
  )
  # Make sure name is correct
  #format.name <- function(x.names) {
  #	gsub("^t\\(x\\)", "", names(x.names))
  #}
  #rownames(temp) <- format.name(rownames(temp))
  temp
}

#' Fits linear model to predict y using each row in x one at a time. Returns coefficients of fit.
#' @param x matrix whose rows are to be regressed vs. y
#' @param y vector of numbers
#' @return slopes of the robust linear model fits of y vs. each variable in x
#' @export
rlm_beta <- function(x, y) {
  if (ncol(x) != length(y)) {
    stop("Must have as many samples in y as in x.")
  }
  # Regress each row vs. y and get the coefficient
  apply(x, 1, function(x.row) {
    tryCatch({
      coef(MASS::rlm(y ~ x.row))[2]
    }, error=function(e) NA)
  })
}

#' Fits linear model to predict y using each row in x and tests whether each slope is significant.
#' @param x matrix whose rows are to be regressed vs. y
#' @param y vector of numbers
#' @param n number of trials for bootstrap and confidence intervals
#' @return slopes of each row in x vs. y, confidence intervals, and P values
#' @examples
#' rlm_beta_test(t(iris[,1:3]), iris[,4], n=100)
#' @importFrom boot boot boot.ci
#' @importFrom stats p.adjust
#' @export
rlm_beta_test <- function(x, y, n=1000) {
  estimate <- rlm_beta(x, y)
  # Select random row indices and do bootstrap
  f <- function(d, ind) {
    rlm_beta(t(d[ind,]), y[ind])
  }
  b.obj <- boot::boot(t(x), f, R=n)
  # Get confidence intervals from bootstraps
  ci <- do.call(rbind, lapply(1:nrow(x), function(i) {
    boot::boot.ci(b.obj, type="perc", index=i)$percent[,4:5]
  }))
  # Permute y and fit models to get P value
  p <- perm_test(function(y) rlm_beta(x, y), y, function() sample(y), n=n)
  fdr <- stats::p.adjust(p, "fdr")
  # Show to final data frame
  data.frame(Estimate = estimate,
             Lower = ci[,1],
             Upper = ci[,2],
             P = p,
             FDR = fdr)
}

#' Returns difference in the regression coefficients of y vs. x in two sample classes
#' @param x matrix of values whose rows are to be regressed vs. y. Columns represent samples.
#' @param y vector of numbers to regress against
#' @param x.class vector of TRUE or FALSE values indicating group assignment of each sample in x.
#' @return vector of difference of regression coefficients for each row in x. Sign convention: TRUE group - FALSE group.
#' @export
rlm_beta_diff <- function(x, y, x.class) {
  if (any(is.na(x.class))) {
    stop("Sample class can't be NA.")
  }
  if (sum(x.class) <= 3 || sum(!x.class) <= 3) {
    return(NA)
  } else {
    beta1 <- rlm_beta(x[,x.class], y[x.class])
    beta2 <- rlm_beta(x[,!x.class], y[!x.class])
    beta1 - beta2
  }
}

#' Test whether there is difference in slopes of y vs. x in two classes of samples
#' @param x matrix where rows = genes, columns = samples
#' @param y value to be regressed with each row
#' @param x.class boolean vector (TRUE or FALSE) for each sample indicating sample group
#' @param n number of permutations
#' @return data frame of estimates, confidence intervals, P values for slopes in each sample class
#'     as well as the difference between the classes.
#' @note In the output, Class1 refers to the TRUE group while Class2 refers to the FALSE group.
#' @importFrom boot boot boot.ci
#' @importFrom stats p.adjust
#' @export
rlm_beta_diff_test <- function(x, y, x.class, n=1000) {
  if (any(is.na(x.class))) {
    stop("Sample class can't be NA.")
  }
  # First calculate for each class
  r1 <- rlm_beta_test(x[,x.class], y[x.class], n=n)
  colnames(r1) <- paste0("Class1_", colnames(r1))
  r2 <- rlm_beta_test(x[,!x.class], y[!x.class], n=n)
  colnames(r2) <- paste0("Class2_", colnames(r2))
  # Difference between the slopes
  est <- rlm_beta_diff(x, y, x.class)
  # Confidence intervals
  f <- function(d, ind) {
    rlm_beta_diff(t(d[ind,]), y[ind], x.class[ind])
  }
  b.obj <- boot::boot(t(x), f, R=n)
  ci <- do.call(rbind, lapply(1:nrow(x), function(i) {
    boot::boot.ci(b.obj, type="perc", index=i)$percent[,4:5]
  }))
  # P value by permuting classes
  p <- perm_test(function(xc) rlm_beta_diff(x, y, xc), x.class, function() sample(x.class), n=n)
  # Output results
  temp <- data.frame(
    Estimate=est,
    Lower=ci[,1],
    Upper=ci[,2],
    P = p,
    FDR = stats::p.adjust(p, "fdr")
  )
  cbind(r1, r2, temp)
}

#' Returns P value of observed value compared to values from null distribution
#' @param observed observed statistic
#' @param null.values a vector of statistics from null distribution
#' @param alternative hypothesis test (can be "two.sided", "less", or "greater")
#' @return p value of observed statistic compared to null distribution
#' @importFrom stats na.omit
#' @examples
#' p_value(1.97, rnorm(1E6))
#' @export
p_value <- function(observed, null.values, alternative="two.sided") {
  if (length(observed) != 1) {
    stop("There must be exactly one observed value for the statistic.")
  }
  if (is.na(observed)) {
    return(NA)
  }
  if (!is.vector(null.values)) {
    stop("Null values must be vector.")
  }
  # Use valid values only
  valid.values <- stats::na.omit(null.values)
  n <- length(valid.values)
  if (alternative == "two.sided") {
    (2 * (min(sum(observed >= null.values, na.rm=TRUE),
              sum(observed <= null.values, na.rm=TRUE))) + 1) / (n + 2)
  } else {
    stop("Not implemented yet. Please use two sided hypothesis only.")
  }
}

#' Tests whether each row in f(x) is equal to 0 using permutation test
#' @param f function that takes x and produces a vector (must have the same names in same order in every trial!)
#' @param x value of x such that f(x) gives observed value
#' @param x.rand function that takes no arguments and gives random value of x under null hypothesis
#' @param n number of trials
#' @return two-sided P value of observed value compared to values assuming null hypothesis
#' @export
perm_test <- function(f, x, x.rand=function() sample(x), n=1000) {
  observed <- f(x)
  null.values <- do.call(cbind, lapply(1:n, function(dummy) {
    f(x.rand())
  }))
  # Return the P values
  sapply(1:length(observed), function(i) {
    p_value(observed[i], null.values[i,])
  })
}

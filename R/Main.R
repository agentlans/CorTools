#' Computes the correlation of each gene with vector of y values
#' @param x matrix of x values (rows = genes, columns = samples)
#' @param y y values (vector with as many samples as x)
#' @param ... extra arguments to pass into correlation function
#' @return matrix of correlations where each row represents correlation of each gene with y
#' @importFrom stats cor
#' @export
cor2 <- function(x, y, ...) {
  if (!is.vector(y)) {
    stop("y must be a vector.")
  }
  if (length(y) != ncol(x)) {
    stop("There must be as many samples in y as in x.")
  }
  # Make x into matrix
  if (is.vector(x)) {
    x <- (matrix(x))
  }
  y <- matrix(y)
  # Compute correlation
  stats::cor(t(x), y, ...)
}
# Example:
#library(multtest)
#data(golub)
#x <- golub[-1,]
#y <- golub[1,]
#cor2(x, y)

#' Returns P value of observed value compared to values from null distribution
#' @param observed observed statistic
#' @param null.values a vector of statistics from null distribution
#' @param alternative hypothesis test (can be "two.sided", "less", or "greater")
#' @return p value of observed statistic compared to null distribution
#' @importFrom stats na.omit
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
    (2 * (min(sum(observed >= null.values),
            sum(observed <= null.values))) + 1) / (n + 2)
  } else {
    stop("Not implemented yet. Please use two sided hypothesis only.")
  }
}
# Example:
# p_value(1.97, rnorm(1E6))

#' Permutation test on correlation between x and y
#' @param x matrix of x values (rows = genes, columns = samples)
#' @param y vector of y values to be correlated with x
#' @param n number of trials for permutation test
#' @param ... extra parameters to pass into cor function
#' @return data frame with components:
#'    Cor (correlation of each row in x with y)
#'    P (P value that correlation is significantly different from 0)
#' @importFrom stats p.adjust
#' @export
cor_test <- function(x, y, n=1000, ...) {
  # Observed value of correlation
  observed <- cor2(x, y, ...)
  # Permute y and get correlation
  null.mat <- do.call(cbind, lapply(1:n, function(dummy) {
    cor2(x, sample(y), ...)
  }))
  # Get the P values
  p.values <- sapply(1:nrow(x), function(i) {
    p_value(observed[i], null.mat[i,])
  })
  # Show confidence intervals and P values
  stat.df <- bootstrap_ci(
    observed, cor_bootstrap(x, y, n=n, ...))
  p.df <- data.frame(
    P = p.values,
    FDR = p.adjust(p.values, "fdr")
  )
  cbind(stat.df, p.df)
}
# Example:
#library(multtest)
#data(golub)
#x <- golub[-1,]
#y <- golub[1,]
#cor_test(x, y)

#' Estimates correlations with y using bootstrapping
#' @param x matrix of x values
#' @param y y values to be correlated against x
#' @param n number of resampling trials
#' @param ... extra arguments to pass to cor
#' @return resampled correlations where each column represents a resampling
#' @export
cor_bootstrap <- function(x, y, n=1000, ...) {
  if (ncol(x) <= 3) {
    warning("Very low number of samples in resampling. Might wish to have more samples or increase resampling fraction.")
  }
  # Do resampling
  resampled.cor <- do.call(cbind, lapply(1:n, function(dummy) {
    # Choose random samples and correlate
    chosen.col <- sample(1:ncol(x), ncol(x), replace=TRUE)
    cor2(x[,chosen.col], y[chosen.col], ...)
  }))
  # Return the resampled correlations
  resampled.cor
}
# library(multtest)
# data(golub)
# x <- golub[-1,]
# y <- golub[1,]
# foo <- cor_bootstrap(x, y)

#' Transform correlation values using Fisher's z transform
#' @param r correlation
#' @return transformed correlation
#' @export
fisher_z <- function(r) {
  0.5 * log((1+r+1E-5)/(1-r+1E-5))
}

#' Returns difference between the Fisher transformed correlations between x and y
#'     in two classes of samples
#' @param x matrix of x values (rows = genes, columns = samples)
#' @param y y values to correlate against (same length as number of samples in x)
#' @param x.classes vector of TRUE and FALSE values for each sample in x, respectively
#' @param ... extra parameters to pass to cor2
#' @return vector of difference between correlations in two sets of samples (TRUE group - FALSE group)
#' @export
two_cor <- function(x, y, x.classes, ...) {
  if (!is.logical(x.classes)) {
    stop("Samples classes must be a vector of TRUE or FALSE values.")
  }
  if (any(is.na(x.classes))) {
    warning("NA in the sample classes may produce unpredictable results.")
  }
  if (length(x.classes) != ncol(x)) {
    stop("Samples classes must be same length as number of samples.")
  }
  # Separate data into to two classes
  x1 <- x[,x.classes]
  y1 <- y[x.classes]
  x2 <- x[,!x.classes]
  y2 <- y[!x.classes]
  # Calculate correlations separately
  # Return difference between Fisher's z transforms
  z1 <- fisher_z(cor2(x1, y1, ...))
  z2 <- fisher_z(cor2(x2, y2, ...))
  z1 - z2
}

#' Bootstrap estimates of difference between two Fisher z-transformed correlations
#' @param x matrix of x values (rows = genes, columns = samples)
#' @param y values to correlate against (must be same length as number of samples in x)
#' @param x.classes vector of TRUE of FALSE values for every sample in x
#' @param n number of trials
#' @param ... parameters to pass to cor2
#' @return matrix of bootstrap estimates of differences in correlations.
#'     Sign convention is TRUE group - FALSE group.
#'     Each column represents bootstrap iteration
#' @export
two_cor_bootstrap <- function(x, y, x.classes, n=1000, ...) {
  do.call(cbind, lapply(1:n, function(dummy) {
    true.col <- which(x.classes == TRUE)
    false.col <- which(x.classes == FALSE)
    # Resample from TRUE and FALSE columns with replacement
    selected <- c(sample(true.col, length(true.col), replace=TRUE),
                  sample(false.col, length(false.col), replace=TRUE))
    # Compute differences in correlations on resampled data frame
    two_cor(x[,selected], y[selected], x.classes[selected], ...)
  }))
}

#' Tests hypothesis that correlation of x with y differs in two sample classes using permutation test
#' @param x matrix of x values (rows = genes, columns = samples)
#' @param y y values to correlate against (must be same length as number of samples in x)
#' @param x.classes a vector of TRUE and FALSE values for every sample in x, respectively
#' @param n number of trials for permutation test
#' @param ... extra arguments to pass into cor
#' @return data frame with components:
#' Cor1 (correlations in samples marked TRUE),
#' Cor2 (correlations in samples marked FALSE)
#' P value that difference in correlations is significant
#' @importFrom stats p.adjust
#' @export
two_cor_test <- function(x, y, x.classes, n=1000, ...) {
  cor.diff <- function(x.classes) two_cor(x, y, x.classes, ...)
  # Get the actual difference in correlations and under null hypothesis
  observed <- cor.diff(x.classes)
  null.values <- do.call(cbind, lapply(1:n, function(dummy) {
    cor.diff(sample(x.classes))
  }))
  # Get P values for every gene
  p.values <- sapply(1:nrow(x), function(i) {
    p_value(observed[i], null.values[i,])
  })
  # Correlation results for C1
  c1.results <- cor_test(x[,x.classes], y[x.classes], ...)
  colnames(c1.results) <- paste0("Cor1_", colnames(c1.results))
  c2.results <- cor_test(x[,!x.classes], y[!x.classes], ...)
  colnames(c2.results) <- paste0("Cor2_", colnames(c2.results))

  # Output correlations of genes in two classes and P value
  cbind(
    c1.results,
    c2.results,
    data.frame(
      DeltaZ = observed,
      P = p.values,
      FDR = p.adjust(p.values, "fdr")
    ))
}

#' Fits y to x using robust linear regression and returns regression coefficients
#' @param x matrix of values where rows = genes, columns = samples
#' @param y value for every sample in x
#' @return coefficients from robust linear regression
#' @importFrom MASS rlm
#' @importFrom stats coef
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
#rlm_coef(t(iris[,1:3]), iris[,4])

#' Tests hypothesis that the coefficients for robust linear regression != 0 using permutation test
#' @param x matrix where genes in rows, samples in columns
#' @param y value for every sample for regression
#' @param n number of trials for permutation test
#' @return data frame containing
#'    Coef coefficients from linear regression
#'    P P values that different from 0
#' @importFrom stats p.adjust
#' @export
rlm_coef_test <- function(x, y, n=1000) {
  # Calculated the observed coefficients and under null hypothesis
  observed <- rlm_coef(x, y)
  null.values <- do.call(cbind, lapply(1:n, function(dummy) {
    rlm_coef(x, sample(y))
  }))
  # P values that coefficients not equal 0 (by permuting samples)
  p <- sapply(1:(nrow(x)+1), function(i) {
    p_value(observed[i], null.values[i,])
  })
  # Get confidence intervals for coefficients (using bootstrap)
  results <- bootstrap_ci(observed, rlm_coef_bootstrap(x, y, n=n))
  cbind(results,
        data.frame(
          P = p,
          FDR = p.adjust(p, "fdr")
        ))
}
# rlm_coef_test(t(iris[,1:3]), iris[,4])

#' Returns robust linear regression coefficients by bootstrapping
#' @param x dataset where rows represent genes, columns represent samples
#' @param y values for each sample in x
#' @param n number of resamples
#' @return matrix of resamples
#' @export
rlm_coef_bootstrap <- function(x, y, n=1000) {
  do.call(cbind, lapply(1:n, function(dummy) {
    # Draw random samples and get coefficients
    chosen <- sample(1:ncol(x), ncol(x), replace=TRUE)
    rlm_coef(x[,chosen], y[chosen])
  }))
}
# Example
#foo <- rlm_coef_bootstrap(t(iris[,1:3]), iris[,4])
#library(reshape)
#boxplot(value~X1, data=melt(foo))
#abline(h=0)

#' Estimates confidence interval based on bootstrap results
#' @param point.est vector of point estimates (each element represents gene)
#' @param bootstrap.resamples matrix of resampled values using bootstrap (each column is a replicate)
#' @param alpha significance level
#' @return data frame containing:
#'     Estimate: point estimates of the statistic
#'     Lower: lower bound of confidence interval
#'     Upper: upper bound of confidence interval
#' @importFrom stats quantile
#' @export
bootstrap_ci <- function(point.est, bootstrap.resamples, alpha=0.05) {
  upper.quantile <- apply(bootstrap.resamples, 1, function(x) quantile(x, 1-alpha/2, na.rm=TRUE))
  lower.quantile <- apply(bootstrap.resamples, 1, function(x) quantile(x, alpha/2, na.rm=TRUE))
  # This is pivot confidence interval
  data.frame(
    Estimate=point.est,
    Lower=2*point.est - upper.quantile,
    Upper=2*point.est - lower.quantile)
}

#' Plots statistics resampled from dataset
#' @param results a matrix where each row represents gene, each column represents resampling trial
#' @param results2 optional second dataset
#' @return ggplot2 graph object
#' @importFrom reshape melt
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_hline
#' @export
plot_resample <- function(results, results2=NULL) {
  #library(ggplot2)
  if (is.null(results2)) {
    # There is only one sample
    graphable <- reshape::melt(results)
    x1 <- factor(graphable$X1)
    value <- graphable$value
    ggplot2::ggplot(graphable, ggplot2::aes(x1, value)) +
      ggplot2::geom_boxplot(width=0.5) +
      ggplot2::geom_hline(ggplot2::aes(yintercept=0), linetype="dashed")
  } else {
    # Melt each dataset separately and record where they came from
    m1 <- reshape::melt(results, id=c("X1", "Dataset"))
    m1$Dataset <- 1
    m2 <- reshape::melt(results2, id=c("X1", "Dataset"))
    m2$Dataset <- 2
    graphable <- rbind(m1, m2)
    # Make plot
    x1 <- factor(graphable$X1)
    value <- graphable$value
    dataset <- factor(graphable$Dataset)
    ggplot2::ggplot(graphable, ggplot2::aes(x1, value)) +
      ggplot2::geom_boxplot(ggplot2::aes(fill=dataset), width=0.5) +
      ggplot2::geom_hline(ggplot2::aes(yintercept=0), linetype="dashed")
  }
}
# g <- plot_resample(rlm_coef_bootstrap(t(iris[,1:3]), iris[,4]))
# g + scale_x_discrete("Variable") + scale_y_continuous("Coefficient")
#
# results <- rlm_coef_bootstrap(t(iris[,1:3]), iris[,4])
# results2 <- rlm_coef_bootstrap(t(iris[,1:3]), iris[,4])
# plot_resample(results, results2)


#' Gets the slope of regression of y vs. x
#' @param x matrix where each row is to be regressed vs. y
#' @param y vector of numbers
#' @importFrom MASS rlm
#' @return Coefficients of regression from each row in x vs. y
#' @export
rlm_beta <- function(x, y) {
  if (ncol(x) != length(y)) {
    stop("Must have as many samples in y as in x.")
  }
  # Regress each row vs. y and get the coefficient
  apply(x, 1, function(x.row) {
    coef(MASS::rlm(y ~ x.row))[2]
  })
}
# rlm_beta(t(iris[,1:3]), iris[,4])

#' Gets bootstrap estimates of slopes of y vs. x
#' @param x matrix where each row is to be regressed vs. y
#' @param y vector of numbers (length same as columns of x)
#' @param n number of trials
#' @return matrix where each row represents slopes for each row in x, each column represents a trial
#' @export
rlm_beta_bootstrap <- function(x, y, n=1000) {
  do.call(cbind, lapply(1:n, function(dummy) {
    # Sample columns with replacement
    ind <- sample(1:ncol(x), ncol(x), replace=TRUE)
    rlm_beta(x[,ind], y[ind])
  }))
}

#' Tests whether slopes of linear regressions equal 0 using permutation test
#' @param x matrix where each row is to be regressed vs. y
#' @param y vector of numbers (length same as columns in x)
#' @param n number of trials
#' @importFrom stats p.adjust
#' @return Data frame of observed value, upper and lower 95% confidence intervals, and P values
rlm_beta_test <- function(x, y, n=1000) {
  observed <- rlm_beta(x, y)
  # Scramble y and calculate slope
  null.values <- do.call(cbind, lapply(1:n, function(dummy) {
    rlm_beta(x, sample(y))
  }))
  # Get P values
  p <- sapply(1:length(observed), function(i) {
    p_value(observed[i], null.values[i,])
  })
  # Bootstrap confidence intervals
  temp <- bootstrap_ci(observed, rlm_beta_bootstrap(x, y, n))
  temp$P <- p
  temp$FDR <- p.adjust(temp$P, "fdr")
  temp
}

#' Tests whether the slopes of y vs. x is the same in two sample groups
#' @param x matrix where each column contains a sample and each row is to be regressed vs. y
#' @param y values to be regressed against
#' @param x.class class of each sample in x, respectively (can be TRUE or FALSE for each sample)
#' @param n number of permutations
#' @return Data frame of bootstrap estimates for each sample group and P value of difference
#'     between two sample groups
#' @importFrom stats p.adjust
#' @export
two_rlm_beta_test <- function(x, y, x.class, n=1000) {
  # Compute difference between slopes in the two classes
  beta_diff <- function(x.class, x.use=x, y.use=y) {
    beta1 <- rlm_beta(x.use[,x.class], y.use[x.class])
    beta2 <- rlm_beta(x.use[,!x.class], y.use[!x.class])
    beta1 - beta2
  }
  # Compute actual difference in slopes
  observed <- beta_diff(x.class)
  null.values <- do.call(cbind, lapply(1:n, function(dummy) {
    beta_diff(sample(x.class))
  }))
  # Get P value of difference
  p <- sapply(1:length(observed), function(i) {
    p_value(observed[i], null.values[i,])
  })
  # Get bootstrapped values for difference
  bootstrapped <- do.call(cbind, lapply(1:n, function(dummy) {
    ind <- sample(1:ncol(x), ncol(x), replace=TRUE)
    beta_diff(x.class[ind], x[,ind], y[ind])
  }))
  beta.diff <- bootstrap_ci(observed, bootstrapped)
  # Do individual tests on the betas
  beta.test1 <- rlm_beta_test(x[,x.class], y[x.class], n)
  beta.test2 <- rlm_beta_test(x[,!x.class], y[!x.class], n)
  colnames(beta.test1) <- paste0("Beta1_", colnames(beta.test1))
  colnames(beta.test2) <- paste0("Beta2_", colnames(beta.test2))
  # Combine together
  temp <- cbind(beta.test1, beta.test2, beta.diff)
  temp$P <- p
  temp$FDR <- p.adjust(p, "fdr")
  temp
}


#' Compares coefficients from two robust linear models using permutation test
#' @param x1 matrix where rows represent predictors and columns represent samples
#' @param x2 exactly like x1 but for a separate model
#' @param y value to predict for each sample
#' @param n number of permutations
#' @return data frame of
#'    Coef1 the coefficients from first model
#'    Coef2 the coefficients from the second model
#'    P P value that coefficients differ between two models
#' @importFrom stats p.adjust
two_rlm_coef_test0 <- function(x1, x2, y, n=1000) {
  # Given y, returns difference between coefficients of models
  # built using x1 and those built using x2
  coef.diff <- function(y) {
    c1 <- rlm_coef(x1, y)
    c2 <- rlm_coef(x2, y)
    common <- intersect(names(c1), names(c2))
    c1[common] - c2[common]
  }
  # Calculate coefficient differences with actual and with y permuted data
  observed <- coef.diff(y)
  null.stats <- do.call(cbind, lapply(1:n, function(dummy) {
    coef.diff(sample(y))
  }))
  p.values <- sapply(1:length(observed), function(i) {
    p_value(observed[i], null.stats[i,])
  })
  # Also do the coefficient tests for each dataset separately
  r1 <- rlm_coef_test(x1, y, n=n)[names(observed)]
  r1$FDR <- p.adjust(r1$P, "fdr") # Readjust P values because fewer variables
  colnames(r1) <- paste0("Model1_", colnames(r1))

  r2 <- rlm_coef_test(x2, y, n=n)[names(observed)]
  r2$FDR <- p.adjust(r2$P, "fdr") # Readjust P values because fewer variables
  colnames(r2) <- paste0("Model2_", colnames(r2))

  # Return coefficients of the two models and P value
  rbind(r1, r2,
        data.frame(
          P = p.values,
          FDR = p.adjust(p.values, "fdr")
        ))
}
#two_rlm_coef_test(t(iris[,1:2]), t(iris[,1:3]), iris[,4])


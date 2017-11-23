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
    stop("Observed value can't be NA.")
  }
  if (!is.vector(null.values)) {
    stop("Null values must be vector.")
  }
  # Use valid values only
  valid.values <- stats::na.omit(null.values)
  n <- length(valid.values)
  if (alternative == "two.sided") {
    2 * (min(sum(observed >= null.values),
            sum(observed <= null.values)) + 1) / (n + 1)
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
  #
  data.frame(
    Cor = observed,
    P = p.values,
    FDR = p.adjust(p.values, "fdr")
  )
}
# Example:
#library(multtest)
#data(golub)
#x <- golub[-1,]
#y <- golub[1,]
#cor_test(x, y)

#' Estimates correlations by taking random subsets of original data and correlating them with y
#' @param x matrix of x values
#' @param y y values to be correlated against x
#' @param frac fraction of samples in x to draw
#' @param n number of resampling trials
#' @param ... extra arguments to pass to cor
#' @return resampled correlations
#' @export
cor_resample <- function(x, y, frac=0.5, n=1000, ...) {
  if (frac*ncol(x) <= 3) {
    warning("Very low number of samples in resampling. Might wish to have more samples or increase resampling fraction.")
  }
  # Do resampling
  resampled.cor <- do.call(cbind, lapply(1:n, function(dummy) {
    # Choose random samples and correlate
    chosen.col <- sample(1:ncol(x), ceiling(frac * ncol(x)))
    cor2(x[,chosen.col], y[chosen.col], ...)
  }))
  # Calculate quantiles on resampled values
  #apply(resampled.cor, 1, function(x) {
  #  quantile(x, c(0.025, 0.975))
  #})
  # Return the resampled correlations
  resampled.cor
}
#library(multtest)
#data(golub)
#x <- golub[-1,]
#y <- golub[1,]
#foo <- cor_resample(x, y)

#' Transform correlation values using Fisher's z transform
#' @param r correlation
#' @return transformed correlation
#' @export
fisher_z <- function(r) {
  0.5 * log((1+r+1E-5)/(1-r+1E-5))
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
  if (!is.logical(x.classes)) {
    stop("Samples classes must be a vector of TRUE or FALSE values.")
  }
  if (any(is.na(x.classes))) {
    warning("NA in the sample classes may produce unpredictable results.")
  }
  if (length(x.classes) != ncol(x)) {
    stop("Samples classes must be same length as number of samples.")
  }
  cor.diff <- function(x.classes) {
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
  # Get the actual difference in correlations and under null hypothesis
  observed <- cor.diff(x.classes)
  null.values <- do.call(cbind, lapply(1:n, function(dummy) {
    cor.diff(sample(x.classes))
  }))
  # Get P values for every gene
  p.values <- sapply(1:nrow(x), function(i) {
    p_value(observed[i], null.values[i,])
  })
  # Output correlations of genes in two classes and P value
  data.frame(
    Cor1 = cor2(x[,x.classes], y[x.classes]),
    Cor2 = cor2(x[,!x.classes], y[!x.classes]),
    P = p.values,
    FDR = p.adjust(p.values, "fdr")
  )
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
  # P values
  p <- sapply(1:(nrow(x)+1), function(i) {
    p_value(observed[i], null.values[i,])
  })
  data.frame(
    Coef = observed,
    P = p,
    FDR = p.adjust(p, "fdr")
  )
}
# rlm_coef_test(t(iris[,1:3]), iris[,4])

#' Returns robust linear regression coefficients by resampling
#' @param x dataset where rows represent genes, columns represent samples
#' @param y values for each sample in x
#' @param frac fraction of samples to use in resample
#' @param n number of resamples
#' @return matrix of resamples
#' @export
rlm_coef_resample <- function(x, y, frac=0.5, n=1000) {
  do.call(cbind, lapply(1:n, function(dummy) {
    # Draw random samples and get coefficients
    chosen <- sample(1:ncol(x), ceiling(frac*ncol(x)))
    rlm_coef(x[,chosen], y[chosen])
  }))
}
# Example
#foo <- rlm_coef_resample(t(iris[,1:3]), iris[,4])
#boxplot(value~X1, data=melt(foo))
#abline(h=0)

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
# g <- plot_resample(rlm_coef_resample(t(iris[,1:3]), iris[,4]))
# g + scale_x_discrete("Variable") + scale_y_continuous("Coefficient")
#
# results <- rlm_coef_resample(t(iris[,1:3]), iris[,4])
# results2 <- rlm_coef_resample(t(iris[,1:3]), iris[,4])
# plot_resample(results, results2)


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
#' @export
two_rlm_coef_test <- function(x1, x2, y, n=1000) {
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
  # Return coefficients of the two models and P value
  data.frame(
    Coef1 = rlm_coef(x1, y)[names(observed)],
    Coef2 = rlm_coef(x2, y)[names(observed)],
    P = p.values,
    FDR = p.adjust(p.values, "fdr")
  )
}
#two_rlm_coef_test(t(iris[,1:2]), t(iris[,1:3]), iris[,4])


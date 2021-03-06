% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{rlm_coef_test}
\alias{rlm_coef_test}
\title{Fits robust linear model on all variables and
    tests whether coefficients of terms in model are equal to 0}
\usage{
rlm_coef_test(x, y, n = 1000)
}
\arguments{
\item{x}{matrix of numbers whose rows are to be regressed with y}

\item{y}{vector of numbers}

\item{n}{number of trials for confidence interval and P values}
}
\value{
coefficients, confidence intervals, and P values from linear model fit
}
\description{
Fits robust linear model on all variables and
    tests whether coefficients of terms in model are equal to 0
}
\examples{
rlm_coef_test(t(iris[,1:3]), iris[,4])
}

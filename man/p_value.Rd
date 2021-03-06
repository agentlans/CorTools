% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{p_value}
\alias{p_value}
\title{Returns P value of observed value compared to values from null distribution}
\usage{
p_value(observed, null.values, alternative = "two.sided")
}
\arguments{
\item{observed}{observed statistic}

\item{null.values}{a vector of statistics from null distribution}

\item{alternative}{hypothesis test (can be "two.sided", "less", or "greater")}
}
\value{
p value of observed statistic compared to null distribution
}
\description{
Returns P value of observed value compared to values from null distribution
}
\examples{
p_value(1.97, rnorm(1E6))
}

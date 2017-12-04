# CorTools

CorTools is an R package for doing permutation tests on correlations and regressions.

Features:
- permutation tests on correlation coefficients and regression coefficients
- estimate correlations and coefficients using resampling
- test whether correlation coefficients are the same in two subsets of samples
- test whether coefficients of variables are the same in different linear models

# Install

```R
# Install prerequisite packages
install.packages(c('MASS', 'boot'))

# Install package from GitHub
library(devtools)
install_github("agentlans/CorTools")
```

# Use
```R
library(CorTools)

# Load example dataset
library(multtest)
data(golub)

# Correlate first 10 genes with the 11th one
cor_test(golub[1:10,], golub[11,])

# Test whether the slope of first 10 genes with 11th one differs in AML vs. ALL
# (To save time, only try 10 random trials)
rlm_beta_diff_test(golub[1:10,], golub[11,], golub.cl == 1, n=10)
```

# License
GPL-3

Software maintainer: Alan Tseng

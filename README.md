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
install.packages(c('ggplot2', 'reshape', 'MASS'))

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

# Plot resampled correlations
library(ggplot2)
plot_resample(cor_bootstrap(golub[1:10,], golub[11,])) +
    scale_x_discrete("Gene") + scale_y_continuous("Correlation with 11th gene")

# Test whether correlations with gene 11 differ in AML vs. ALL
two_cor_test(golub[1:10,], golub[11,], golub.cl == 1)

# Plot correlations in the two classes
r1 <- cor_bootstrap(golub[1:10, golub.cl==1], golub[11, golub.cl==1])
r2 <- cor_bootstrap(golub[1:10, golub.cl==0], golub[11, golub.cl==0])
plot_resample(r1, r2) + scale_x_discrete("Gene") +
    scale_y_continuous("Correlation with 11th gene")
```

# License
GPL-3

Software maintainer: Alan Tseng

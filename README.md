
<!-- README.md is generated from README.Rmd. Please edit that file -->

# H2x2Factorial

<!-- badges: start -->
<!-- badges: end -->

This `H2x2Factorial` package implements the sample size methods for
hierarchical 2x2 factorial trials with unequal cluster sizes. The sample
size calculations support five types of hypothesis tests: (A1) test for
marginal cluster-level treatment effect, (A2) test for marginal
individual-level treatment effect, (B) interaction test for the two
treatments, (C) joint test for the two marginal treatment effects, (D)
intersection union test for the two marginal treatment effects.
Finite-sample considerations are included for the tests involving the
marginal cluster-level treatment effect, due to the degree of freedom
issues. Three functions are currently contained for predicting the power
or sample size based on given design parameters as well as delivering
illustrative tables or line plots. Specifically, the
`calc.H2x2Factorial` function calculates required number of clusters for
a specific test to achieve a given power, or predicts the actual power
given specified sample size resources, with or without finite-sample
considerations. The `table.H2x2Factorial` function creates a data frame
to show a series of sample size predictions by providing varying mean
cluster sizes, intraclass correlation coefficients, or coefficient of
variations of cluster sizes (CV). The `graph.H2x2Factorial` function
plots sample size requirements under different CV in the form of the
combinations of mean cluster sizes and number of clusters. All of the
hypothesis tests and sample size methodologies are formalized in “Sample
Size Calculation in Hierarchical 2x2 Factorial Trials with Unequal
Cluster Sizes” (under review).

## Installation

The development version of the packages can be installed from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("BillyTian/H2x2Factorial")
```

## Example

Here shows some basic examples:

``` r
library(H2x2Factorial)

calc.H2x2Factorial(n.input=10, delta_x=0.2, delta_z=0.1, test="joint", correction=T, seed.mix=123456, CV=0.38, rho=0.1)

table.H2x2Factorial(delta_x=0.2, delta_z=0.1, m_bar=c(10,50,100), CV=c(0, 0.3, 0.5), rho=c(0.01, 0.1), test="cluster")

graph.H2x2Factorial(power=0.9, test="cluster", rho=0.1)
```

# robustAR
This R package provides the implementation of the robust AR statistic introduced in Klooster and Zhelonkin (202X).

## Installation
To install `robustAR`, use the function `install_github()` of the `devtools` package.

```R
install.packages("devtools")
devtools::install_github("Jens-Klooster/robustAR", build_vignettes = TRUE)
library(robustAR)
```

After the installation, the robustAR vignette provides further examples of how to use the function `robustAR`. The vignette can be opened as follows

```R
browseVignettes(package="robustAR")
```
and by then clicking on HTML.

## Basic Example

```R
# True beta value
beta <- 1

# amount of observations we generate
n <- 250

# we generate correlated error terms
eps1 <- rnorm(n)
eps2 <- 0.5*eps1 + 0.5*rnorm(n)

# we generate a control variable
W <- rnorm(n)

# we generate an instrumental variable
Z <- rnorm(n) 

# first-stage equation
X <- 0.15 * Z + 1 * W + eps1

# second-stage equation
Y <- beta * X + 2 * W + eps2
```

The robust AR statistic based on Tukey's biweight function and hat matrix weights can then be computed as follows
```R
robustAR(Y = Y, X = X, Z = Z, W = W, beta0 = 1, 
         type = "Tukey", 
         weighting = "hat")
```

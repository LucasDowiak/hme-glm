# hme-glm
An R implementation of the hierarchical mixture-of-experts model proposed by
Jordan and Jacobs (1994).

## Getting Started
This package is not yet available through CRAN. It needs to be cloned from the 
Github repo and built from source. The code below is for a generic linux shell but
it assumes you have installed R and it is accessible via the command line.

```bash
git clone https://github.com/LucasDowiak/hme-glm.git
cd hme-glm
R CMD build .
R CMD install hmeglm_0.0.0.1.tar.gz
```

Note that the name of the binary file (e.g. 'hmeglm_0.0.0.1.tar.gz') will change
for every new update and release so please check before running `R CMD install`.
Barring any encountered errors, the package is now available and ready to load in
an R sessions.

Invoking the model is simple and adheres to the R language's normal syntax. For
model specification, there are two important inputs. The input variable `tree`
controls the architecture of the gating network and the regression specification
is governed by a `formula` object. The general structure of the `formula`is
'y ~ x | z' where 'y' is the dependent variable, 'x' the set of explanatory
variables in the local experts regressions and 'z' is the set of variables in the
gating network that are used to allocate (in a soft-max manner) input patterns to
the most appropriate local regression. Let's demonstrate this mixture model on
the Iris data set:

```r
library(hmeglm)
data(iris)
tree_ <- c("0", "0.1", "0.2") # Specify the tree architecture of the gating network
mod <- hme(tree=tree_,
           formula="Sepal.Width ~ Petal.Width | Petal.Width",
           data=iris,
           family=gaussian())
summary(mod)
```

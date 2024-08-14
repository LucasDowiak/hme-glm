# hme-glm
An R implementation of the hierarchical mixture-of-experts (HME) model proposed by
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

## Model API

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
mod <- hme(tree=c("0", "0.1", "0.2"),                         # Set the architecture of the gating network
           formula="Sepal.Width ~ Petal.Width | Petal.Width", # Specify the regression
           data=iris,
           family=gaussian())
summary(mod)
```

## Outline of the Theory

The HME is presented as a standard mixture model. For a given input and
output pair $(x_t, y_t)$, each expert provides a probabilistic
model relating input row $\textbf{x}_{t}$ to output $y_{t}$:

$$\begin{equation}
  P^{m}_{t} \equiv P^{m}(y_{t}|\textbf{x}_{t}; \boldsymbol{\beta}^{m}), \quad m = 1,2,...,M
\end{equation}$$

where $m$ is one of the $M$ component experts in the mixture. The experts are
combined with associated weights into a mixture distribution:

$$
\begin{equation}
  P(y_{t} | \textbf{x}_{t}; \, \boldsymbol{\beta}) = \sum_{m=1}^{M} \Pi(m|t) P^{m}(y_{t} | \textbf{x}_{t}; \boldsymbol{\beta}^{m})
\end{equation}
$$

Here, $\Pi(m|t)$ is the probability that the input unit $t$ belongs to expert $m$
and has the usual restrictions: $0 \leq \Pi(m|t) \leq 1$ for each $m$ and
$\sum_{m} \Pi(m|t) = 1$. The novel aspect of the HME is the functional form it
uses to model $\Pi(m|t)$. It uses a set of "gates", structured as a tree, to
recursively partition the input space and apply a set of local regression to the
partitioned space. To do so, the gating network includes a second set of co-variates $\textbf{z}_{t}$ and parameter vector $\boldsymbol{\Omega}$:

$$
\begin{equation}
  P(y_{t} | \textbf{x}_{t}, \textbf{z}_{t}; \, \boldsymbol{\beta}, \boldsymbol{\Omega}) = \sum_{m=1}^{M} \Pi(m | \textbf{z}_{t}; \boldsymbol{\Omega}) P^{m}(y_{t} | \textbf{x}_{t}; \boldsymbol{\beta}^{m})
\end{equation}
$$

The graphic below demonstrates various ways to structure the gating network so that
it results in four expert mixture. Gating nodes are recorded as circles while
expert nodes are denoted by squares. The labels of each node are critical to running
the model and are collectively passed to the argument `tree` as a character vector.

![The gating architectures **A** and **B** are "symmetric" trees. While **A** uses
a single multinomial split to allocate inputs, architecture **B** uses a set of
recursive binary splits. Architectures **C** and **D** are both asymmetric trees but
differ in their application of asymmetry. Gating network **C** uses recursive binary
splits that branch from the left side of the tree. Network **D** first uses a tertiary
split at the first level and then a single binary split at the second level.](./images/gating_architectures.png)


# diffindiff: Nonparametric Heterogeneous Treatment Effect Estimation in Repeated Cross Sectional Designs

This package implements methods for non-parametric difference-in-differences estimation, as proposed by Nie, Lu, and Wager (2021) as an update version to the previous paper ``Robust Nonparametric Difference-in-Differences Estimation" (Lu, Nie, and Wager 2019). We consider a
setup where we observe data `(X, Y, S, T)` generated according
to the following model:
```
X ~ P(X)
S ~ P(S|X), T ~ P(T|X), where S, T takes values {0,1}
Y = b(X) + T*f(X) + S*h(X) + T*S*tau(X) + noise
```
and we want to estimate the average effect `TAU = E[tau(X)]`. The general framework is flexible: auxilary parameters such as `b(X)`, `f(X)` and `h(x)` can be estimated with any preferred machine learning method. Our methods are implemented with the lasso (glmnet), and `cv.glmnet` is used to perform cross-fitting.

### Installation

To install this package in R, run the following commands:
```R
library(devtools) 
install_github("xnie/diffindiff")
```
### Example usage

```R
library(diffindiff)
library(glmnet)
library(rlearner)
library(EQL)
library(magrittr)

n = 100; p=12
X = matrix(rnorm(n * p), n, p)
tau = 1
f = X[,3]
h = X[,4]
S.treat.prob = 0.4
T.treat.prob = 0.4
Si = rbinom(n, 1, S.treat.prob)
Ti = rbinom(n, 1, T.treat.prob)
b = pmax(X[,1] + X[,2], 0)
Y = b + Ti*f + Si*h + Ti*Si*tau + rnorm(n)

# generate a basis
make_matrix = function(x) stats::model.matrix(~.-1, x)
X = data.frame(X) %>% make_matrix
order = 3
Basis = generate.basis(X,order)

# If we assume that the treatment effect tau is constant, which is true
# for the setup above, then we can run the transformed regression
TR_fit = DiD(X = Basis,
  Y = Y,
  Ti = Ti,
  Si = Si,
  constant_eff = "constant")
  
# If we do not assume tau is constant, we can run DiD-AMLE
gamma.minimax = minimax(Basis, Ti, Si)
DiD_AMLE_fit = DiD(X = Basis,
  Y = Y,
  Ti = Ti,
  Si = Si,
  constant_eff = "non_constant",
  gamma = gamma.minimax)
  
# To repeat the heterogeneous treatment effect simulation experiments 
cd experiments
bash simulations/start_hte.sh
```

#### References
Xinkun Nie, Chen Lu, Stefan Wager.
<b>Nonparametric Heterogeneous Treatment Effect Estimation in Repeated Cross Sectional Designs</b>
2021.
[<a href="https://arxiv.org/abs/1905.11622v2">arxiv</a>]

# CVXQUAD

CVXQUAD is a collection of functions to be used with the MATLAB-based convex optimization tool [CVX](http://www.cvxr.com/cvx/). It implements a new approximation strategy to treat the exponential cone as well as various functions based on matrix logarithm using symmetric cone solvers. This package is based on the paper:

```
Semidefinite approximations of matrix logarithm
Hamza Fawzi, James Saunderson and Pablo A. Parrilo.
```

# Installation

Download the zip file https://github.com/hfawzi/cvxquad/archive/master.zip and add it to your MATLAB path.

## Replacing successive approximation
To replace the successive approximation functionality of CVX whenever the exponential cone is used (e.g., when using rel_entr or in GP mode), copy the file "exponential/exponential.m" to the folder "sets" in your CVX installation (you may want to keep a copy of the existing file in case you want to revert to the successive approximation method).

# Example

The following code uses the ```quantum_rel_entr``` function of CVXQUAD to compute the nearest correlation matrix to a given matrix M, in the quantum relative entropy sense.

```
n = 4;
M = randn(n,n);
M = M*M';
cvx_begin
  variable X(n,n) symmetric
  minimize quantum_rel_entr(M,X)
  subject to
    diag(X) == ones(n,1)
cvx_end
```

# Functions and sets

| Function | | |
| --- | --- | --- |
| rel_entr_quad(x,y) | x.*log(x./y) | convex in (x,y) |
| quantum_entr(X) | -trace(X*logm(X)) | concave in X |
| quantum_rel_entr(X,Y) | trace(X*(logm(X)-logm(Y))) | convex in (X,Y) |
| trace_logm(X,C) | trace(C*logm(X)) | concave in X<br />(C fixed positive semidefinite matrix) |
| trace_mpower(X,t,C) | trace(C*X^t) |  concave in X for t in [0,1]<br /> convex in X for t in [-1,0] or [1,2]<br />(C fixed positive semidefinite matrix) |
| lieb_ando(X,Y,K,t) | trace(K' \* X^{1-t} \* K \* Y^t) |  concave in (X,Y) for t in [0,1]<br /> convex in (X,Y) for t in [-1,0] or [1,2]<br /> (K is a fixed matrix)|

| Set | |
| --- | --- |
| op_rel_entr_epi_cone | Operator relative entropy cone |
| matrix_geo_mean_hypo_cone | Matrix geometric mean hypograph cone |
| matrix_geo_mean_epi_cone | Matrix geometric mean epigraph cone |

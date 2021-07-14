# CVXQUAD

CVXQUAD is a collection of functions to be used with the MATLAB-based convex optimization tool [CVX](http://www.cvxr.com/cvx/). It augments CVX with various convex/concave functions based on matrix logarithm such as the von Neumann entropy, or the quantum relative entropy (see below for list of functions). This package is based on the paper:

```
Semidefinite approximations of matrix logarithm
Hamza Fawzi, James Saunderson and Pablo A. Parrilo
```

available at https://arxiv.org/abs/1705.00812.

# News

*July 14, 2021*: Added support for approximations that are formal upper or lower bounds on log -- the Pade approximations from the paper above are neither upper nor lower bounds. The new approximations are based on Gauss-Radau quadrature, see [doc/log_approx_bounds.pdf](doc/log_approx_bounds.pdf) for more details. Functions quantum_entr, quantum_rel_entr, trace_logm and op_rel_entr_epi_cone now have an additional parameter, 'apx', which can be used to specify which approximation to be used (apx=-1,0,+1).

# Installation

Unpack the zip file https://github.com/hfawzi/cvxquad/archive/master.zip and add the folder to your MATLAB path.

## Replacing successive approximation

*Update*: CVX 2.2 can now handle the exponential cone natively when used with Mosek 9. See [CVX's documentation](http://web.cvxr.com/cvx/doc/advanced.html#the-successive-approximation-method). If you have access to Mosek, no need to use any approximations.

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

# Citing

To cite the package in your work, you can use the following bibtex code:

```
@article{cvxquad,
  title={Semidefinite approximations of the matrix logarithm},
  author={Fawzi, Hamza and Saunderson, James and Parrilo, Pablo A.},
  year={2018},
  journal={Foundations of Computational Mathematics},
  note={Package cvxquad at \url{https://github.com/hfawzi/cvxquad}}
}
```

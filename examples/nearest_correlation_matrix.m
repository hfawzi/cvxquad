addpath('../');

% Nearest correlation matrix in the quantum relative entropy sense

n = 4;
M = randn(n,n);
M = M*M';
cvx_begin
  variable X(n,n) symmetric
  minimize quantum_rel_entr(M,X)
  diag(X) == ones(n,1)
cvx_end
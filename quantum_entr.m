function cvx_optval = quantum_entr(X,m,k)

%QUANTUM_ENTR    Quantum (von Neumann) entropy
%   QUANTUM_ENTR(X) returns -trace(X*logm(X)) where the logarithm is base e
%   (and not base 2!). X must be a positive semidefinite matrix.
%
%   Disciplined convex programming information:
%      QUANTUM_ENTR is concave in X.
%      This function implements the semidefinite programming approximation
%      given in the reference below. Parameters m and k control the
%      accuracy of this approximation: m is the number of quadrature nodes
%      to use and k the number of square-roots to take. See reference for
%      more details. Default (m,k) = (3,3).
%
%   REQUIRES: op_rel_entr_epi_cone
%   Implementation uses H(X) = -trace( D_{op}(X||I) ) where D_{op} is the
%   operator relative entropy:
%       D_{op}(X||Y) = X^{1/2}*logm(X^{1/2} Y^{-1} X^{1/2})*X^{1/2}
%
%AUTHORS
%   Hamza Fawzi, James Saunderson and Pablo A. Parrilo
%
%REFERENCE
%   This code is based on the paper: "Semidefinite approximations of matrix
%   logarithm" by Hamza Fawzi, James Saunderson and Pablo A. Parrilo

% Check input arguments
if nargin == 1
    m = 3;
    k = 3;
end
if size(X,1) ~= size(X,2)
    error('Input must be a square matrix.');
end

if isnumeric(X)
    cvx_optval = -quantum_rel_entr(X,eye(size(X,1)));
elseif cvx_isconstant(X)
    cvx_optval = cvx(quantum_entr(cvx_constant(X),m,k));
elseif cvx_isaffine(X)
    n = size(X,1);
    iscplx = ~isreal(X);
    cvx_begin
        if iscplx
            variable TAU(n,n) hermitian
        else
            variable TAU(n,n) symmetric
        end
        {X,eye(n),TAU} == op_rel_entr_epi_cone(n,iscplx,m,k);
        maximize (-trace(TAU));
    cvx_end
else
    error('Disciplined convex programming error:\n The input has to be an affine expression');
end

end
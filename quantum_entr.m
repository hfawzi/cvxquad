function cvx_optval = quantum_entr(X,m,k,apx)

%QUANTUM_ENTR    Quantum (von Neumann) entropy
%   QUANTUM_ENTR(X) returns -trace(X*logm(X)) where the logarithm is base e
%   (and not base 2!). X must be a positive semidefinite matrix.
%
%   Disciplined convex programming information:
%      QUANTUM_ENTR is concave in X.
%      This function implements the semidefinite programming approximation
%      given in the reference below.
%
%      Parameters m and k control the accuracy of this approximation:
%      m is the number of quadrature nodes to use and k the number of
%      square-roots to take. See reference for more details.
%      Default (m,k) = (3,3).
%
%      Parameter apx indicates which approximation r of the von Neumann
%      entropy to use:
%      - apx = +1: Upper approximation of entropy (H(X) <= r(X))
%      - apx = -1: Lower approximation of entropy (r(X) <= H(X))
%      - apx = 0 (Default): Pade approximation (neither upper nor lower),
%                           but slightly better accuracy than apx=+1 or -1.
%      The upper and lower approximation are based on rational functions
%      derived from Gauss-Radau quadrature, see documentation in the 'doc'
%      folder.
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

if nargin < 4
    % By default use Pade approximant
    apx = 0;
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
        % Since H(X) = -tr[D_{op}(X||I)], to get an upper/lower bound (resp.) on
        % H(X), we need a lower/upper bound (resp.) on D_{op}. That's why
        % we use -apx and not apx
        {X,eye(n),TAU} == op_rel_entr_epi_cone(n,iscplx,m,k,eye(n),-apx);
        maximize (-trace(TAU));
    cvx_end
else
    error('Disciplined convex programming error:\n The input has to be an affine expression');
end

end
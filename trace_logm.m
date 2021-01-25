function cvx_optval = trace_logm(X,C,m,k)

%TRACE_LOGM    Trace of logarithm
%   TRACE_LOGM(X) returns trace(logm(X)) where X is a positive definite
%   matrix.
%   TRACE_LOGM(X,C) returns trace(C*logm(X)) where C is a positive
%   semidefinite matrix of the same size as X.
%
%   Disciplined convex programming information:
%      TRACE_LOGM(X,C) is concave in X, provided C is a (fixed) positive
%      semidefinite matrix.
%      This function implements the semidefinite programming approximation
%      given in the reference below.  Parameters m and k control the
%      accuracy of this approximation: m is the number of quadrature nodes
%      to use and k the number of square-roots to take. See reference for
%      more details. Default (m,k) = (3,3).
%
%   REQUIRES: op_rel_entr_epi_cone
%   Implementation uses trace(C*logm(X)) = -trace(C*D_{op}(I||X))
%   where D_{op} is the operator relative entropy:
%     D_{op}(X||Y) = X^{1/2}*logm(X^{1/2} Y^{-1} X^{1/2})*X^{1/2}
%
%AUTHORS
%   Hamza Fawzi, James Saunderson and Pablo A. Parrilo
%
%REFERENCE
%   This code is based on the paper: "Semidefinite approximations of matrix
%   logarithm" by Hamza Fawzi, James Saunderson and Pablo A. Parrilo

if nargin < 1 || nargin > 4
    error('Wrong number of arguments');
end

if ~ismatrix(X) || size(X,1) ~= size(X,2)
    error('X must be a square matrix');
end
if nargin < 2 || numel(C) == 0
    C = eye(size(X,1));
else
    if size(C,1) ~= size(C,2) || size(C,1) ~= size(X,1)
        error('C must be a positive semidefinite matrix of the same size as X');
    end
    C = (C+C')/2;
    e = eig(C);
    tol = 1e-9;
    if min(e) < -tol
        error('C must be positive semidefinite');
    end
end
if nargin == 2
    m = 3;
    k = 3;
end

if isnumeric(X)
    cvx_optval = -quantum_rel_entr(C,X)+quantum_rel_entr(C,eye(size(C,1)));
elseif cvx_isconstant(X)
    cvx_optval = cvx(trace_logm(cvx_constant(X),C));
elseif cvx_isaffine(X)
    n = size(X,1);
    iscplx = ~isreal(X) || ~isreal(C);
    cvx_begin
        if iscplx
            variable TAU(n,n) hermitian
        else
            variable TAU(n,n) symmetric
        end
        % -logm(X) <= TAU
        {eye(n),X,TAU} == op_rel_entr_epi_cone(n,iscplx,m,k);
        maximize -trace(C*TAU)
    cvx_end
else
    error('Disciplined convex programming error:\n The input has to be an affine expression');
end

end

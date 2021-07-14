function cvx_optval = quantum_rel_entr(A,B,m,k,apx)

%QUANTUM_REL_ENTR    Quantum relative entropy
%   QUANTUM_REL_ENTR(A,B) returns trace(A*(logm(A)-logm(B))) where A and B
%   are positive semidefinite matrices such that \im(A) \subseteq \im(B)
%   (otherwise the function evaluates to infinity). Note this function uses
%   logarithm base e (and not base 2!).
%
%   Disciplined convex programming information:
%      QUANTUM_REL_ENTR(A,B) is convex (jointly) in (A,B)
%      This function implements the semidefinite programming approximation
%      given in the reference below.
%
%      Parameters m and k control the accuracy of this approximation:
%      m is the number of quadrature nodes to use and k the number of
%      square-roots to take. See reference for more details.
%      Default (m,k) = (3,3).
%
%      Parameter apx indicates which approximation r of the relative
%      entropy function to use:
%      - apx = +1: Upper approximation (D(A|B) <= r(A,B))
%      - apx = -1: Lower approximation (r(A,B) <= D(A|B))
%      - apx = 0 (Default): Pade approximation (neither upper nor lower),
%                           but slightly better accuracy than apx=+1 or -1.
%      The upper and lower approximation are based on rational functions
%      derived from Gauss-Radau quadrature, see documentation in the 'doc'
%      folder.
%
%   REQUIRES: op_rel_entr_epi_cone
%   Implementation uses the expression
%      D(A||B) = e'*D_{op} (A \otimes I || I \otimes B) )*e
%   where D_{op} is the operator relative entropy and e = In(:) where In is
%   the nxn identity matrix.
%
%AUTHORS
%   Hamza Fawzi, James Saunderson and Pablo A. Parrilo
%
%REFERENCE
%   This code is based on the paper: "Semidefinite approximations of matrix
%   logarithm" by Hamza Fawzi, James Saunderson and Pablo A. Parrilo

if nargin < 2
    error('Not enough input arguments');
end

if ~ismatrix(A) || ~ismatrix(B) || size(A,1) ~= size(A,2) ...
        || size(B,1) ~= size(B,2) || size(A,1) ~= size(B,1)
    error('A and B must be square matrices of the same size');
end
if nargin == 2
    m = 3;
    k = 3;
end

if nargin < 5
    % By default use Pade approximant
    apx = 0;
end

if isnumeric(A) && isnumeric(B)
    % Compute trace(A*(logm(A)-logm(B)))
    tol = 1e-9;
    % Check that A and B are symmetric/Hermitian
    if norm(A-A','fro') > tol*norm(A,'fro') ...
            || norm(B-B','fro') > tol*norm(B,'fro')
        error('A and B must be symmetric or Hermitian matrices');
    end
    % Enforce A and B are symmetric/Hermitian
    A = (A+A')/2; B = (B+B')/2;
    % Diagonalize A and B and check that im(A) \subseteq im(B)
    [V,a] = eig(A); a = diag(a);
    [W,b] = eig(B); b = diag(b);
    if min(a) < -tol || min(b) < -tol
        error('A and B must be positive semidefinite');
    end
    ia = find(a > tol);
    ib = find(b > tol);
    u = a'*abs(V'*W).^2;
    if any(u(b <= tol) > tol)
        error('D(A||B) is infinity because im(A) is not contained in im(B)');
    else
        r1 = sum(a(ia).*log(a(ia))); % trace(A*logm(A))
        r2 = u(ib)*log(b(ib)); % trace(A*logm(B))
        cvx_optval = r1 - r2;
    end
elseif cvx_isconstant(A)
    cvx_optval = -quantum_entr(A,m,k) - trace_logm(B,A,m,k,-apx);
elseif cvx_isconstant(B)
    cvx_optval = -quantum_entr(A,m,k,-apx) - trace(A*logm(B));
elseif cvx_isaffine(A) && cvx_isaffine(B)
    n = size(A,1);
    In = eye(n);
    e = In(:);
    iscplx = ~isreal(A) || ~isreal(B);
    cvx_begin
        variable tau;
        {kron(A,eye(n)),kron(eye(n),conj(B)),tau} == op_rel_entr_epi_cone(n^2,iscplx,m,k,e,apx);
        minimize tau;
    cvx_end
else
    error('Disciplined convex programming error:\n The input has to be an affine expression');
end

end
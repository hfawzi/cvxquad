function cvx_optval = trace_mpower(A,t,C)

%TRACE_MPOWER    Trace of matrix power
%   trace_mpower(A,t) returns trace(A^t) where A is a positive semidefinite
%   matrix and t \in [-1,2].
%   trace_mpower(A,t,C) returns trace(C*A^t) where C is a positive
%   semidefinite matrix.
%
%   Disciplined convex programming information:
%      When t \in [0,1], TRACE_MPOWER(A,t,C) is concave in A (for fixed
%      positive semidefinite matrix C) and convex for t \in [-1,0] or
%      [1,2].
%
%AUTHORS
%   Hamza Fawzi and James Saunderson
%
%REFERENCE
%   This code is based on the paper: "Lieb's concavity theorem, matrix
%   geometric means and semidefinite optimization" by Hamza Fawzi and James
%   Saunderson (arXiv:1512.03401)

if nargin < 2
    error('Not enough input arguments');
end
if ~ismatrix(A) || size(A,1) ~= size(A,2)
    error('X must be a square matrix');
end
if nargin < 3
    C = eye(size(A,1));
else
    % Check that C is positive semidefinite
    C = (C+C')/2;
    if any(eig(C) < -1e-6)
        error('C has to be positive semidefinite');
    end
end

if isnumeric(A)
    cvx_optval = trace(C*mpower(A,t));
elseif cvx_isaffine(A)
    if t < -1 || t > 2
        error('t must be between -1 and 2');
    end
    n = size(A,1);
    iscplx = ~isreal(A);
    cvx_begin
        if iscplx
            variable T(n,n) hermitian
        else
            variable T(n,n) symmetric
        end
        if t >= 0 && t <= 1
            % Concave function
            maximize trace(C*T)
            {eye(n),A,T} == matrix_geo_mean_hypo_cone(n,t,iscplx,0);
        else
            % Convex function
            minimize trace(C*T)
            {eye(n),A,T} == matrix_geo_mean_epi_cone(n,t,iscplx,0);
        end
    cvx_end    
else
    error('Disciplined convex programming error:\n The input has to be an affine expression');
end

end
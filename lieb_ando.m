function cvx_optval = lieb_ando(A,B,K,t)

%LIEB_ANDO    Lieb's function
%   LIEB_ANDO(A,B,K,t) returns  trace(K' * A^{1-t} * K * B^t) where A and B
%   are positive semidefinite matrices and K is an arbitrary matrix
%   (possibly rectangular).
%
%   Disciplined convex programming information:
%      LIEB_ANDO(A,B,K,t) is concave in (A,B) for t \in [0,1], and convex
%      in (A,B) for t in [-1,0] or [1,2]. K is a fixed matrix.
%
%AUTHORS
%   Hamza Fawzi and James Saunderson
%
%REFERENCE
%   This code is based on the paper: "Lieb's concavity theorem, matrix
%   geometric means and semidefinite optimization" by Hamza Fawzi and James
%   Saunderson (arXiv:1512.03401)

if nargin ~= 4
    error('Wrong number of arguments');
end

if ~ismatrix(A) || size(A,1) ~= size(A,2)
    error('A must be a square matrix');
end

if ~ismatrix(B) || size(B,1) ~= size(B,2)
    error('B must be a square matrix');
end

if size(K,1) ~= size(A,1) || size(K,2) ~= size(B,1)
    error('K must have the same number of rows as A and the same number of columns as B');
end

if isnumeric(A) && isnumeric(B)
    cvx_optval = real(trace(K'*mpower(A,1-t)*K*mpower(B,t)));
elseif isnumeric(A)
    KAK = K'*mpower(A,1-t)*K;
    KAK = (KAK+KAK')/2;
    cvx_optval = trace_mpower(B,t,KAK);
elseif isnumeric(B)
    KBK = K*mpower(B,t)*K';
    KBK = (KBK+KBK')/2;
    cvx_optval = trace_mpower(A,1-t,KBK);
elseif cvx_isaffine(A) && cvx_isaffine(B)
    n = size(A,1);
    m = size(B,1);
    iscplx = ~isreal(A) || ~isreal(B) || ~isreal(K);
    Kvec = reshape(K.',n*m,1);
    KvKv = Kvec*Kvec';
    KvKv = (KvKv+KvKv')/2;
    cvx_begin
        if iscplx
            variable T(n*m,n*m) hermitian
        else
            variable T(n*m,n*m) symmetric
        end
        if t >= 0 && t <= 1
            % Concave function
            maximize (trace(KvKv*T))
            {kron(A,eye(m)),kron(eye(n),conj(B)),T} == matrix_geo_mean_hypo_cone(n*m,t,iscplx,0);
            %{A,conj(B),T} == lieb_hypo_cone(n,m,t,iscplx);
        elseif (t >= -1 && t <= 0) || (t >= 1 && t <= 2)
            % Convex function
            minimize (trace(KvKv*T))
            {kron(A,eye(m)),kron(eye(n),conj(B)),T} == matrix_geo_mean_epi_cone(n*m,t,iscplx,0);
            %{A,conj(B),T} == ando_epi_cone(n,m,t,iscplx);
        else
            error('t must be between -1 and 2');
        end
    cvx_end
else
    error('Disciplined convex programming error:\n The input has to be an affine expression');
end

end
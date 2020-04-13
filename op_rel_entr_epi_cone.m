function cvx_optpnt = op_rel_entr_epi_cone(sz,iscplx,m,k,e)

%OP_REL_ENTR_EPI_CONE    Operator relative entropy cone
%   OP_REL_ENTR_EPI_CONE(sz,iscplx) returns a CVX triple {X,Y,TAU} of
%   matrices constrained to satisfy
%       X^{1/2} * logm(X^{1/2}*Y^{-1}*X^{1/2}) * X^{1/2} <= TAU
%   where:
%     * The inequality "<=" is in the positive semidefinite order
%     * X,Y,TAU are symmetric (Hermitian if iscplx=1) matrices of size sz
%     * logm is the matrix logarithm (in base e)
%
%   Usage example:
%     variable X(4,4) symmetric
%     variable Y(4,4) symmetric
%     variable TAU(4,4) symmetric
%     {X,Y,TAU} == op_rel_entr_epi_cone(4,0)
%
%   Note on parameter sz:
%     It is possible to provide for sz an array [sz(1) sz(2) sz(3) ...].
%     This will return an array of prod(sz(3:end)) triples, where each
%     triple has size sz(1)=sz(2) and is constrained to live in the
%     operator relative entropy ocne. Note that sz(1) and sz(2) must be
%     equal in this case. See documentation for CVX's "sets/semidefinite.m"
%     for more information.
%
%   OP_REL_ENTR_EPI_CONE(sz,iscplx,m,k)  The additional parameters m and k
%   can be used to control the accuracy of the approximation (see reference
%   for more details). m is the number of quadrature nodes to use and k is
%   the number of square-roots to take. Default (m,k) = (3,3).
%
%   OP_REL_ENTR_EPI_CONE(sz,iscplx,m,k,e)  If an additional matrix e is
%   provided of size nxr (where n=sz(1)), the returned tuple(s) {X,Y,TAU}
%   is constrained to satisfy instead:
%      e'*(D_{op}(X||Y))*e <= TAU.
%   Note that TAU here is of size rxr. The default case corresponds to
%   e=eye(n). When r is small compared to n this can be helpful to reduce
%   the size of small LMIs from 2nx2n to (n+r)x(n+r).
%
%AUTHORS
%   Hamza Fawzi, James Saunderson and Pablo A. Parrilo
%
%REFERENCE
%   This code is based on the paper: "Semidefinite approximations of matrix
%   logarithm" by Hamza Fawzi, James Saunderson and Pablo A. Parrilo


if nargin < 1 || nargin > 5
    error('Wrong number of arguments');
end

if nargin < 2
    iscplx = 0;
end
if nargin == 2
    m = 3;
    k = 3;
end


if length(sz) == 1
    sz = [sz sz];
end
if size(sz,2) == 1
    sz = sz';
end
if sz(1) ~= sz(2)
    error('The first two elements of sz must be equal');
end

if nargin < 5
    e = eye(sz(1));
else
    if size(e,1) ~= sz(1)
        error('The number of rows of e must be the same as sz(1)');
    end
end
r = size(e,2);

% Compute Gauss-Legendre quadrature nodes and weights on [0,1]
[s,w] = glquad(m);

% Make sure that w is a column vector
if size(w,1) == 1
    w = w';
end

cvx_begin set
    if iscplx
        variable X(sz) hermitian
        variable Y(sz) hermitian
        variable Z(sz) hermitian
        variable TAU([r r sz(3:end)]) hermitian
        variable T([r r sz(3:end) m]) hermitian
    else
        variable X(sz) symmetric
        variable Y(sz) symmetric
        variable Z(sz) symmetric
        variable TAU([r r sz(3:end)]) symmetric
        variable T([r r sz(3:end) m]) hermitian
    end
    
    % X #_{1/2^k} Y >= Z
    {X,Y,Z} == matrix_geo_mean_hypo_cone(sz,1/(2^k),iscplx,0);

    for ii=1:m
        sr.type = '()';
        sr.subs = repmat({':'},1,length(sz));
        sr.subs{length(sz)+1} = ii;
        % subsref(T,sr) corresponds to T(:,...,:,ii)
        % Also note that we are dividing by w here because it is easier
        % to do this than to do sum w_i T(:,...,:,ii) later (cf. line that
        % involves TAU)
        
        [e'*X*e - s(ii)*subsref(T,sr)/w(ii)   e'*X;
            X*e       (1-s(ii))*X+s(ii)*Z] == semidefinite([r+sz(1) r+sz(1) sz(3:end)],iscplx);
        
    end
    (2^k)*sum(T,length(sz)+1) + TAU == semidefinite([r r sz(3:end)],iscplx);
cvx_end

cvx_optpnt = cvxtuple(struct('X',X,'Y',Y,'TAU',TAU));

end

function [s,w] = glquad(m)
    % Compute Gauss-Legendre quadrature nodes and weights on [0,1]
    % Code below is from [Trefethen, "Is Gauss quadrature better than
    % Clenshaw-Curtis?", SIAM Review 2008] and computes the weights and
    % nodes on [-1,1].
    beta = .5./sqrt(1-(2*(1:(m-1))).^(-2)); % 3-term recurrence coeffs
    T = diag(beta,1) + diag(beta,-1); % Jacobi matrix
    [V,D] = eig(T); % eigenvalue decomposition
    s = diag(D); % nodes
    w = 2*V(1,:).^2; % weights
    % Translate and scale to [0,1]
    s = (s+1)/2; w = w'/2;
end
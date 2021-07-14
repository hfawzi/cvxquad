function cvx_optpnt = op_rel_entr_epi_cone(sz,iscplx,m,k,e,apx)

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
%   OP_REL_ENTR_EPI_CONE(sz,iscplx,m,k)  The additional parameters m
%   and k can be used to control the accuracy of the approximation (see
%   reference for more details). m is the number of quadrature nodes to use
%   and k is the number of square-roots to take. Default (m,k) = (3,3).
%
%   OP_REL_ENTR_EPI_CONE(sz,iscplx,m,k,e)  If an additional matrix e is
%   provided of size nxr (where n=sz(1)), the returned tuple(s) {X,Y,TAU}
%   is constrained to satisfy instead:
%      e'*(D_{op}(X||Y))*e <= TAU.
%   Note that TAU here is of size rxr. The default case corresponds to
%   e=eye(n). When r is small compared to n this can be helpful to reduce
%   the size of small LMIs from 2nx2n to (n+r)x(n+r).
%
%   OP_REL_ENTR_EPI_CONE(sz,iscplx,m,k,e,apx)  The last parameter apx
%   indicates which approximation of the logarithm to use.
%   Possible values:
%   - apx = +1: Upper approximation on D_{op} [inner approximation of the
%               operator relative entropy cone]
%   - apx = -1: Lower approximation on D_{op} [outer approximation of the
%               operator relative entropy cone]
%   - apx = 0 (Default): Pade approximation (neither upper nor lower) but
%                        slightly better accuracy than apx = +1 or -1.
%
%AUTHORS
%   Hamza Fawzi, James Saunderson and Pablo A. Parrilo
%
%REFERENCE
%   This code is based on the paper: "Semidefinite approximations of matrix
%   logarithm" by Hamza Fawzi, James Saunderson and Pablo A. Parrilo


if nargin < 1 || nargin > 6
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

if nargin < 6
    % By default use Pade approximant
    apx = 0;
end


% Note that D_{op}(x,y) is the (operator) perspective of -log, namely
% D_{op}(x,y) = x*(-log)(y/x).
% If f <= log <= F then we get upper and lower bound on D_{op} as follows:
%   -x*F(y/x) <= D_{op}(x,y) <= -x*f(y/x)
% This means:
%   upper bound on D_{op} <-> lower bound on log
%   lower bound on D_{op} <-> upper bound on log

% Compute Gauss-Legendre quadrature nodes and weights on [0,1]
if apx == +1
    % Upper bound on D_{op} <-> lower bound on log
    % Use Gauss-Radau quadrature with an endpoint at 1.
    % Note: m=1 gives the lower bound 1-1/x <= log(x).
    [s,w] = gaussradau(m,1);
elseif apx == -1
    % Lower bound on D_{op} <-> upper bound on log
    % Use Gauss-Radau quadrature with an endpoint at 0. This leads to 
    % the (m,m-1) Pade approximant of log.
    % Note: m=1 gives the upper bound log(x) <= x-1.
    [s,w] = gaussradau(m,0);
elseif apx == 0
    % Pade approximant: use Gaussian quadrature
    [s,w] = glquad(m);
else
    error('Unknown value of parameter apx. (Should be either +1, -1, or 0)');
end

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
        
        if s(ii) < 1e-6
            % for s=0, x*f_s(z/x) = z - x
            e'*(Z-X)*e - subsref(T,sr)/w(ii) == semidefinite([r r sz(3:end)],iscplx);
        else
            % x f_s(z/x) >= T_i / w_i
            [e'*X*e - s(ii)*subsref(T,sr)/w(ii)   e'*X;
                X*e       (1-s(ii))*X+s(ii)*Z] == semidefinite([r+sz(1) r+sz(1) sz(3:end)],iscplx);
        end
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

function [s,w] = gaussradau(m,endpoint)
    % Computes m-point Gauss-Radau quadrature nodes and weights on [0,1]
    % endpoint should be either 0 or 1
    % We use the modified Golub-Welsch method, explained in the first two
    % pages of
    % https://www.cs.purdue.edu/homes/wxg/selected_works/section_08/164.pdf
    % The nodes are computed as the eigenvalues of the modified Jacobi
    % matrix T, where T(m,m) is set to be
    %   T(m,m) = a - beta(m)*P_{m-1}(a) / P_m(a)
    % where P_m are the monic orthogonal polynomial with respect to the
    % base measure on [-1,1] (here uniform measure).
    % In our case P_m = Legendre polynomials which satisfies
    %       P_m(1) = 2^m / (2m choose m)
    % so that P_m(1) / P_{m-1}(1) = (2m-3)/(m-1).
    if m == 1
        s = endpoint;
        w = 1;
    else
        beta = .5./sqrt(1-(2*(1:(m-1))).^(-2));
        a = sign(endpoint-.5)*(1 - beta(m-1)^2*(2*m-3)/(m-1));
        T = diag(beta,1) + diag(beta,-1); % Jacobi matrix
        T(m,m) = a;
        [V,D] = eig(T); % eigenvalue decomposition
        s = diag(D);
        w = 2*V(1,:).^2;
        % Translate and scale to [0,1]
        s = (s+1)/2; w = w'/2;
    end
    % The following must be true:
    % * One of the s_i's is equal to endpoint
    % * Furthermore:
    %       sum_{i=1}^m w_i s_i^j = 1/(j+1) for all j=0,...,2*m-2
    %     (  i.e., w'*s.^(0:(2*m-2)) == 1./(1+(0:(2*m-2)))  )
end
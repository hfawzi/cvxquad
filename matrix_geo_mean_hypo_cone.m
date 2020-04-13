function cvx_optpnt = matrix_geo_mean_hypo_cone(sz,t,iscplx,fullhyp)

%MATRIX_GEO_MEAN_HYPO_CONE    Matrix geometric mean cone.
%   MATRIX_GEO_MEAN_HYPO_CONE(sz,t,iscplx,1) returns a CVX tuple {A,B,T} of
%   matrices constrained to satisfy
%      A #_{t} B >= T
%   where:
%     * A #_{t} B is the t-weighted geometric mean of A and B (see
%     reference below for details). Parameter t should be in [0,1].
%     * The inequality ">=" is in the positive semidefinite order
%     * A,B,T are symmetric (Hermitian if iscplx=1) matrices of size sz
%
%   Note on parameter sz:
%     It is possible to provide for sz an array [sz(1) sz(2) sz(3) ...].
%     This will return an array of prod(sz(3:end)) triples, where each
%     triple has size sz(1)=sz(2) and is constrained to live in the
%     operator relative entropy ocne. Note that sz(1) and sz(2) must be
%     equal in this case. See documentation for CVX's "sets/semidefinite.m"
%     for more information.
%
%   Note on parameter fullhyp:
%     In many applications one doesn't need the full hypograph
%       hyp_t = {(A,B,T) : A #_{t} B >= T}
%     but rather it is enough to work with a convex set C_t that satisfies
%       (A,B,A #_{t} B) \in C_t
%       (A,B,T) \in C_t  =>  A #_{t} B >= T
%     In this case one should set fullhyp = 0. The SDP description will be
%     (slightly) smaller. (By default fullhyp is set to 1).
%
%AUTHORS
%   Hamza Fawzi and James Saunderson
%
%REFERENCE
%   This code is based on the paper: "Lieb's concavity theorem, matrix
%   geometric means and semidefinite optimization" by Hamza Fawzi and James
%   Saunderson (arXiv:1512.03401)

if nargin < 2 || nargin > 4
    error('Wrong number of arguments');
end

if nargin < 2
    t = 1/2;
end

if nargin < 3
    iscplx = 0;
end

if nargin < 4
    fullhyp = 1;
end

if t < 0 || t > 1
    error('t has to be between 0 and 1');
end

if length(sz) == 1
    sz = [sz sz];
end
if size(sz,2) == 1
    sz = sz';
end
dsz = [2*sz(1) 2*sz(2) sz(3:end)];

[p,q] = rat(t);
isPowerOfTwo = @(k) ~bitand(k,k-1);

cvx_begin set
    if iscplx
        variable A(sz) hermitian
        variable B(sz) hermitian
        variable T(sz) hermitian
    else
        variable A(sz) symmetric
        variable B(sz) symmetric
        variable T(sz) symmetric
    end
    
    
    % Reduce to the case fullhyp=0 first
    if fullhyp
        if iscplx
            variable W(sz) hermitian
        else
            variable W(sz) symmetric
        end
        {A,B,W} == matrix_geo_mean_hypo_cone(sz,t,iscplx,0);
        W - T == semidefinite(sz,iscplx);
    else
        % Now we assume fullhyp=0, the recursion are simpler in this case
        if t == 0
            A == semidefinite(sz,iscplx);
            B == semidefinite(sz,iscplx);
            A - T == 0; % fullhyp = 0
        elseif t == 1
            A == semidefinite(sz,iscplx);
            B == semidefinite(sz,iscplx);
            B - T == 0; % fullhyp = 0
        elseif t == 1/2
            [A T; T B] == semidefinite(dsz,iscplx);
        elseif isPowerOfTwo(q)
            % Dyadic number
            if iscplx
                variable Z(sz) hermitian
            else
                variable Z(sz) symmetric
            end
            if t < 1/2
                [A T; T Z] == semidefinite(dsz,iscplx);
                {A,B,Z} == matrix_geo_mean_hypo_cone(sz,2*t,iscplx,0);
            else
                [B T; T Z] == semidefinite(dsz,iscplx);
                {A,B,Z} == matrix_geo_mean_hypo_cone(sz,2*t-1,iscplx,0);
            end
        elseif isPowerOfTwo(p) && t > 1/2
            % Numerator is a power of two and t > 1/2
            if iscplx
                variable Z(sz) hermitian
            else
                variable Z(sz) symmetric
            end
            {A,T,Z} == matrix_geo_mean_hypo_cone(sz,(2*p-q)/p,iscplx,0);
            [Z T; T B] == semidefinite(dsz,iscplx);
        elseif t < 1/2
            % General case t < 1/2
            % Decompose t in t = (p/2^l) * (2^l/q) where l=floor(log2(q))
            l = floor(log2(q));
            if iscplx
                variable X(sz) hermitian
            else
                variable X(sz) symmetric
            end
            {A,B,X} == matrix_geo_mean_hypo_cone(sz,p/(2^l),iscplx,0);
            {A,X,T} == matrix_geo_mean_hypo_cone(sz,(2^l)/q,iscplx,0);
        else
            % General case t \in [1/2,1]
            % Use transformation t <-> 1-t to bring it to case t \in [0,1/2]
            {B,A,T} == matrix_geo_mean_hypo_cone(sz,1-t,iscplx,0);
        end
    end
cvx_end

cvx_optpnt = cvxtuple(struct('A',A,'B',B,'T',T));

end
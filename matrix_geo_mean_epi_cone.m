function cvx_optpnt = matrix_geo_mean_epi_cone(sz,t,iscplx,fullhyp)

%MATRIX_GEO_MEAN_EPI_CONE    Matrix geometric mean cone.
%   MATRIX_GEO_MEAN_EPI_CONE(sz,t,iscplx) returns a CVX tuple {A,B,T} of
%   matrices constrained to satisfy
%      A #_{t} B <= T
%   where:
%     * A #_{t} B is the t-weighted geometric mean of A and B (see
%     reference below for details). Parameter t should be in [-1,0] or
%     [1,2].
%     * The inequality "<=" is in the positive semidefinite order
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
%   Note on fullhyp: matrix_geo_mean_epi_cone will always return a full
%   epigraph cone (unlike matrix_geo_mean_hypo_cone) and so this parameter
%   is not really used. (It is here just for consistency with the hypo_cone
%   function.)
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

if t < -1 || (t > 0 && t < 1) || t > 2
    error('t has to be in [-1,0] or [1,2]');
end

if length(sz) == 1
    sz = [sz sz];
end
if size(sz,2) == 1
    sz = sz';
end
dsz = [2*sz(1) 2*sz(2) sz(3:end)];

cvx_begin set
    if iscplx
        variable A(sz) hermitian
        variable B(sz) hermitian
        variable T(sz) hermitian
        variable Z(sz) hermitian
    else
        variable A(sz) symmetric
        variable B(sz) symmetric
        variable T(sz) symmetric
        variable Z(sz) symmetric
    end
    if t <= 0
        [T A; A Z] == semidefinite(dsz,iscplx);
        {A,B,Z} == matrix_geo_mean_hypo_cone(sz,-t,iscplx,0);
    elseif t >= 1
        [T B; B Z] == semidefinite(dsz,iscplx);
        {A,B,Z} == matrix_geo_mean_hypo_cone(sz,2-t,iscplx,0);
    end
cvx_end

cvx_optpnt = cvxtuple(struct('A',A,'B',B,'T',T));

end
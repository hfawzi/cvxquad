function cvx_optval = rel_entr_quad(x,y,m,k)
%REL_ENTR_QUAD Relative entropy
%   REL_ENTR_QUAD(x,y) returns x.*log(x./y) where x and y are a vector of
%   positive numbers
%
%   Disciplined convex programming information:
%      REL_ENTR_QUAD(x,y) is convex in (x,y)
%      This function implements a second-order cone approximation
%      given in the reference below.  Parameters m and k control the
%      accuracy of this approximation: m is the number of quadrature nodes
%      to use and k the number of square-roots to take. See reference for
%      more details. Default (m,k) = (3,3).
%
%   REQUIRES: op_rel_entr_epi_cone
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

% -- This portion is based on cvx's rel_entr.m --
% x and y are allowed to be vectors.
% Determine size of output sz.
sx = size( x ); xs = all( sx == 1 );
sy = size( y ); ys = all( sy == 1 );
if xs,
    sz = sy;
    x = repmat(x,sz);
elseif ys,
    sz = sx;
    y = repmat(y,sz);
elseif isequal( sx, sy ),
    sz = sx;
else
    error( 'Dimensions of x and y are not compatible.' );
end
% ----------------------------------------------

if nargin < 4
    m = 3;
    k = 3;
end

if isnumeric(x) && isnumeric(y)
    cvx_optval = rel_entr(x,y);
elseif cvx_isconstant(x) && cvx_isconstant(y)
    cvx_optval = cvx(rel_entr_quad(cvx_constant(X),m,k));
elseif cvx_isaffine(x) || cvx_isaffine(y)
    cvx_begin separable
        variable z(sz)
        minimize z
        {shiftdim(x,-2),shiftdim(y,-2),shiftdim(z,-2)} == ...
            op_rel_entr_epi_cone([1 1 sz],0,m,k);
    cvx_end
end

end


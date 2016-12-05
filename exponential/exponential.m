function cvx_optpnt = exponential( sx )

%EXPONENTIAL   Exponential cone.
%   Implementation of the exponential cone based on approximations given in
%   the reference below. This file takes the same arguments as CVX's
%   original exponential.m file.
%
%AUTHORS
%   Hamza Fawzi, James Saunderson and Pablo A. Parrilo
%   (based on the original exponential.m file in CVX)
%
%REFERENCE
%   This code is based on the paper: "Semidefinite approximations of matrix
%   logarithm" by Hamza Fawzi, James Saunderson and Pablo A. Parrilo


error( nargchk( 0, 1, nargin ) ); %#ok

%
% Check size vector
%

if nargin == 0 || isempty( sx ),
    sx = [1,1]; %#ok
else
    [ temp, sx ] = cvx_check_dimlist( sx, true ); %#ok
    if ~temp,
        error( 'First argument must be a dimension vector.' );
    end
end


%
% Build the cone
%

cvx_begin set
    variables x( sx ) y( sx ) z( sx )
    m = 3; % Pade level
    k = 3; % number of square roots
    fprintf('=====================================\n');
    fprintf('Using Pade approximation for exponential\n');
    fprintf('cone with parameters m=%d, k=%d\n',m,k);
    fprintf('=====================================\n');
    {shiftdim(y,-2),shiftdim(z,-2),shiftdim(-x,-2)} == ...
        op_rel_entr_epi_cone([1 1 sx],0,m,k);
cvx_end

cvx_optpnt = cvxtuple( struct( 'x', x, 'y', y, 'z', z ) );

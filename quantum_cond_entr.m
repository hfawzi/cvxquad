function cvx_optval = quantum_cond_entr(rho,dim,sys)

%QUANTUM_COND_ENTR    Quantum conditional entropy
%   If rho is a symmetric (or Hermitian) matrix of size na*nb then
%     quantum_cond_entr(rho,[na nb]) = H(A|B)
%
%   Use the third parameter sys to get H(B|A) instead:
%     quantum_cond_entr(rho,[na nb],2) = H(B|A)
%
%   Concave function of rho.
% 
%   NOTE: This function requires the quantinf package by Toby Cubitt (for
%   the partial trace), see http://www.dr-qubit.org/Matlab_code.html

if ~exist('TrX','file')
    error('This function requires the quantinf package (http://www.dr-qubit.org/Matlab_code.html)\n');
end

if nargin == 2
    sys = 1;
end

if sys == 1
    cvx_optval = -quantum_rel_entr(rho,kron(eye(dim(1)),TrX(rho,sys,dim)));
elseif sys == 2
    cvx_optval = -quantum_rel_entr(rho,kron(TrX(rho,sys,dim),eye(dim(2))));
end

end
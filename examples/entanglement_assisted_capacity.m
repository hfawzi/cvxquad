% Entanglement-assisted classical capacity of a quantum channel

% Dimensions of input, output, and environment spaces of channel
na = 2; nb = 2; ne = 2;
% AD(gamma) = isometry representation of amplitude damping channel
AD = @(gamma) [1 0; 0 sqrt(gamma); 0 sqrt(1-gamma); 0 0];
U = AD(0.2);
assert(size(U,1) == nb*ne && size(U,2) == na);

cvx_begin sdp
    variable rho(na,na) hermitian;
    maximize ((quantum_cond_entr(U*rho*U',[nb ne]) + ...
                    quantum_entr(TrX(U*rho*U',2,[nb ne])))/log(2));
    subject to
        rho >= 0; trace(rho) == 1;
cvx_end

% Quantum capacity of degradable channels

% Example: amplitude damping channel
% na = channel input dimension
% nb = channel output dimension
% ne = channel environment dimension
% nf = degrading map environment dimension
na = 2; nb = 2; ne = 2; nf = 2;

% AD(gamma) = isometry representation of amplitude damping channel
AD = @(gamma) [1 0; 0 sqrt(gamma); 0 sqrt(1-gamma); 0 0];
gamma = 0.2;
U = AD(gamma);

% Unitary representation of degrading map
W = AD ((1-2*gamma)/(1-gamma));

% Ic(rho) = coherent information (see Eq. (13) of paper)
Ic = @(rho) quantum_cond_entr( ...
                W*applychan(U,rho,'isom',[na nb])*W', [ne nf], 2)/log(2);

% Quantum capacity = maximum of Ic (for degradable channels)
cvx_begin sdp
    variable rho(na,na) hermitian
    maximize (Ic(rho));
    rho >= 0; trace(rho) == 1;
cvx_end
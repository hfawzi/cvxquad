% Compute capacity of a cq-channel

% Example 2.16 in 
%   "Efficient Approximation of Quantum Channel Capacities" by Sutter et al. (arXiv:1407.8202)
rho1 = [1 0; 0 0]; H1 = quantum_entr(rho1);
rho2 = (1/2) * [1 1; 1 1]; H2 = quantum_entr(rho2);

cvx_begin
    variables p1 p2;
    maximize ((quantum_entr(p1*rho1 + p2*rho2) - p1*H1 - p2*H2)/log(2))
    subject to
        p1 >= 0; p2 >= 0; p1+p2 == 1;
cvx_end

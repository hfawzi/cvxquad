% Compute lower bound on relative entropy of entanglement (PPT relaxation)

na = 2; nb = 2;
rho = randRho(na*nb); % Generate a random bipartite state rho

cvx_begin sdp
    variable tau(na*nb,na*nb) hermitian;
    minimize (quantum_rel_entr(rho,tau)/log(2));
    tau >= 0; trace(tau) == 1;
    Tx(tau,2,[na nb]) >= 0; % Positive partial transpose constraint
cvx_end

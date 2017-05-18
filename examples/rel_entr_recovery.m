% Compute relative entropy of recovery of a tripartite state rhoABC,
% defined as:
%    min_{L : B->BC}    D(rhoABC || (id_A \otimes L)(rhoAB))
% where the minimization is over trace-preserving completely positive maps
% L : B->BC

if ~exist('TrX','file') || ~exist('sysexchange','file') ...
        || ~exist('tensor','file') || ~exist('applychan','file')
    error('This code requires the quantinf package (http://www.dr-qubit.org/Matlab_code.html)\n');
end

% Random state
na = 2; nb = 2; nc = 2;
rhoABC = randRho(na*nb*nc);

% An explicit state where I(A:C|B) < rel. entr. of recovery
% na = 2; nb = 2; nc = 2;
% theta = 0.3;
% psiAC1 = [1; 0; 0; 0];
% psiAC2 = [0; cos(theta); sin(theta); 0];
% psiABC = sysexchange((kron([1;0],psiAC1) + kron([0;1],psiAC2))/sqrt(2), [1 2], [nb na nc]);
% rhoABC = psiABC*psiABC';
% assert(abs(trace(rhoABC)-1) <= 1e-9);

rhoAB = TrX(rhoABC, 3, [na nb nc]);

% Compute relative entropy of recovery
cvx_begin sdp
    % Channel \Lambda:B->BC (represented by its Choi-Jamiolkowski matrix)
    variable L_B_BC(nb^2*nc,nb^2*nc) hermitian
    
    % chanout_ABC will hold (id_A \otimes \Lambda)(rhoAB)
    variable chanout_ABC(na*nb*nc, na*nb*nc) hermitian
    
    % Objective function
    minimize (quantum_rel_entr(rhoABC, chanout_ABC) / log(2));
    
    % \Lambda must be trace-preserving and completely positive
    L_B_BC >= 0;
    TrX(L_B_BC,2,[nb nb*nc]) == eye(nb);

    % Form Choi matrix of id_A \otimes \Lambda given that of \Lambda
    % w is maximally entangled state on A
    w = zeros(na^2,1); w([1:na+1:na^2]) = 1; w = w*w';
    L_AB_ABC = sysexchange(tensor(w,L_B_BC),[2 3],[na na nb nb*nc]);

    % chanout_ABC is the result of applying id_A\otimes\Lambda to rhoAB
    chanout_ABC == applychan(L_AB_ABC,rhoAB,'choi2',[na*nb na*nb*nc]);
    
cvx_end

% Conditional mutual information
rhoBC = TrX(rhoABC,1,[na nb nc]); assert(size(rhoBC,1) == nb*nc);
rhoAB = TrX(rhoABC,3,[na nb nc]); assert(size(rhoAB,1) == na*nb);
rhoB = TrX(rhoABC,[1 3],[na nb nc]); assert(size(rhoB,1) == nb);
IACcondB = (quantum_entr(rhoBC) + quantum_entr(rhoAB) - quantum_entr(rhoABC) - quantum_entr(rhoB))/log(2);

fprintf('I(A:C|B) = %.4f\n', IACcondB);
fprintf('Rel. entr. of recovery = %.4f\n', cvx_optval);

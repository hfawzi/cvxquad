% Test symmetry reduction for constraints of the form
%   P_g (A \otimes I, I \otimes B) \psd T
% where g is operator concave

n = 4;

% Symmetric trace-less matrices of size n
Psym1 = [];
Pskewsym1 = [];
for ii=1:n
    for jj=(ii+1):n
        Eij = zeros(n,n);
        Eij(ii,jj) = 1;
        Eij(jj,ii) = 1;
        Psym1 = [Psym1 Eij(:)];
        Fij = zeros(n,n);
        Fij(ii,jj) = 1;
        Fij(jj,ii) = -1;
        Pskewsym1 = [Pskewsym1 Fij(:)];
    end
end
for ii=2:n
    EE = zeros(n,n);
    EE(1,1) = 1;
    EE(ii,ii) = -1;
    Psym1 = [Psym1 EE(:)];
end
In = eye(n);
Pw = In(:);
Pw = Pw/norm(Pw);

nsym = n*(n+1)/2-1;
nskewsym = n*(n-1)/2;

assert(size(Psym1,1) == n^2 && size(Psym1,2) == nsym);
assert(size(Pskewsym1,1) == n^2 && size(Pskewsym1,2) == nskewsym);

[Psym,~] = qr(Psym1,0);
[Pskewsym,~] = qr(Pskewsym1,0);

% 

A = randn(n,n); A = A*A';
B = randn(n,n); B = B*B';

C = randn(n^2,n^2); C = C*C';

cvx_begin

    variable T(n^2,n^2) symmetric
    
    {Pw'*kron(A,In)*Pw, Pw'*kron(In,B)*Pw, Pw'*T*Pw} == matrix_geo_mean_hypo_cone(1,1/2,0);
    {Psym'*kron(A,In)*Psym, Psym'*kron(In,B)*Psym, Psym'*T*Psym} == matrix_geo_mean_hypo_cone(nsym,1/2,0);
    {Pskewsym'*kron(A,In)*Pskewsym, Pskewsym'*kron(In,B)*Pskewsym, Pskewsym'*T*Pskewsym} == matrix_geo_mean_hypo_cone(nskewsym,1/2,0);
    
    Pw'*T*Psym == 0;
    Pw'*T*Pskewsym == 0;
    Psym'*T*Pskewsym == 0;
    
    maximize (Pw'*T*Pw);
    
cvx_end

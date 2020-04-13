function cvxquad_tests()

    % CVXQUAD_TESTS
    % Run some unit tests on cvxquad's functions

    rng(0);

    test_matrix_geo_mean_hypo_cone();
    test_op_rel_entr_epi_cone();
    test_quantum_rel_entr();
    test_lieb_ando();
    
end

function test_matrix_geo_mean_hypo_cone()
    fprintf('---- TESTING matrix_geo_mean_hypo_cone ----\n');
    nvec = [3 5];
    tvec = [1/2 1/4 1/8 1/16 3/4 7/8 15/16 2/3 6/7];
    for n=nvec
        for t=tvec
            for cplx=[0 1]
                A = randPSD(n,cplx);
                B = randPSD(n,cplx);
                cvx_quiet(true);
                cvx_begin
                    if cplx
                        variable T(n,n) hermitian
                    else
                        variable T(n,n) symmetric
                    end
                    maximize trace(T)
                    {A,B,T} == matrix_geo_mean_hypo_cone(n,t,cplx,0);
                cvx_end
                ABt = A^(1/2)*mpower(A^(-1/2)*B*A^(-1/2),t)*A^(1/2);
                [p,q] = rat(t);
                fprintf('n=%d, t=%d/%d, cplx=%d: ', n,p,q,cplx);
                assert(norm(ABt-T) <= 1e-6,'Test failed matrix_geo_mean_hypo_cone n=%d, t=%d/%d, cplx=%d', n,p,q,cplx);
                fprintf('OK\n');
            end
        end
    end
    
    % Additional test for fullhyp=1
    A = [6.25 0; 0 16];
    B = [2 1; 1 2];
    % These matrices satisfy A^(1/2) >= B and so (A,eye(2),B) \in hyp_{1/2}
    % A and B however do not satisfy A >= B^2 and so [A B; B eye(2)] not
    % psd, and using fullhyp=0 will give an infeasible SDP.
    cvx_quiet(true)
    cvx_begin
      minimize 0
      {A,eye(2),B} == matrix_geo_mean_hypo_cone(2,1/2,0,1);
    cvx_end
    fprintf('Testing fullhyp=1: ')
    assert(cvx_optval == 0, 'Test matrix_geo_mean_hypo_cone with fullhyp=1 failed');
    fprintf('OK\n');
    
    
    fprintf('matrix_geo_mean_hypo_cone test OK\n');
end

function test_op_rel_entr_epi_cone()
    fprintf('---- TESTING op_rel_entr_epi_cone ----\n');
    nvec = [3 5 10];
    for n=nvec
        for cplx=[0 1]
            A = randPSD(n,cplx); A = A/trace(A);
            B = randPSD(n,cplx); B = B/trace(B);
            cvx_quiet(true);
            cvx_begin
                if cplx
                    variable T(n,n) hermitian
                else
                    variable T(n,n) symmetric
                end
                minimize trace(T)
                {A,B,T} == op_rel_entr_epi_cone(n,cplx);
            cvx_end
            DopAB = A^(1/2)*logm(A^(1/2)*inv(B)*A^(1/2))*A^(1/2);
            fprintf('n=%d, cplx=%d: ',n,cplx);
            % The 1e-3 was set by inspection
            assert(norm(DopAB-T) <= 1e-3,'Test failed op_rel_entr_epi_cone n=%d, cplx=%d, error=%.4e',n,cplx,norm(DopAB-T));
            fprintf('OK\n');
        end
    end
    
    fprintf('op_rel_entr_epi_cone test OK\n');
end

function test_quantum_rel_entr()
    fprintf('---- TESTING quantum_rel_entr ----\n');
    nvec = [2 3];
    for n=nvec
        for cplx=[0 1]
            A = randPSD(n,cplx); A = A/trace(A);
            B = randPSD(n,cplx); B = B/trace(B);
            
            cvx_clear;
            cvx_quiet(true);
            cvx_begin
                if cplx
                    variable X(n,n) hermitian
                else
                    variable X(n,n) symmetric
                end
                minimize (quantum_rel_entr(X,B))
                X == A;
            cvx_end
            val1 = cvx_optval;
            
            cvx_clear
            cvx_quiet(true);
            cvx_begin
                if cplx
                    variable Y(n,n) hermitian
                else
                    variable Y(n,n) symmetric
                end
                minimize (quantum_rel_entr(A,Y))
                Y == B;
            cvx_end
            val2 = cvx_optval;
            
            cvx_clear
            cvx_quiet(true);
            cvx_begin
                if cplx
                    variable X(n,n) hermitian
                    variable Y(n,n) hermitian
                else
                    variable X(n,n) symmetric
                    variable Y(n,n) symmetric
                end
                minimize (quantum_rel_entr(X,Y))
                X == A;
                Y == B;
            cvx_end
            val12 = cvx_optval;
            
            DAB = quantum_rel_entr(A,B);
            fprintf('n=%d, cplx=%d: ',n,cplx);
            
            % The 1e-3 was set by inspection here
            assert(abs(DAB-val1) <= 1e-3,'Test failed quantum_rel_entr n=%d, cplx=%d, error=%.4e',n,cplx,abs(DAB-val1));
            assert(abs(DAB-val2) <= 1e-3,'Test failed quantum_rel_entr n=%d, cplx=%d, error=%.4e',n,cplx,abs(DAB-val2));
            assert(abs(DAB-val12) <= 1e-3,'Test failed quantum_rel_entr n=%d, cplx=%d, error=%.4e',n,cplx,abs(DAB-val12));
            fprintf('OK\n');
        end
    end
    
    fprintf('quantum_rel_entr test OK\n');
end

function test_lieb_ando()

    fprintf('---- TESTING lieb_ando ----\n');
    
    nvec = [2 3];
    tvec = [1/2 1/4 3/4 1/8 3/2 5/4];
    for n=nvec
        for t=tvec
            for cplx=[0 1]
                In = eye(n);
                A = randPSD(n,cplx); A = A/trace(A);
                B = randPSD(n,cplx); B = B/trace(B);

                cvx_clear;
                cvx_quiet(true);
                cvx_begin
                    if cplx
                        variable X(n,n) hermitian
                    else
                        variable X(n,n) symmetric
                    end
                    if t >= 0 && t <= 1
                        maximize (lieb_ando(X,B,In,t))
                    else
                        minimize (lieb_ando(X,B,In,t))
                    end
                    X == A;
                cvx_end
                val1 = cvx_optval;

                cvx_clear
                cvx_quiet(true);
                cvx_begin
                    if cplx
                        variable Y(n,n) hermitian
                    else
                        variable Y(n,n) symmetric
                    end
                    if t >= 0 && t <= 1
                        maximize (lieb_ando(A,Y,In,t))
                    else
                        minimize (lieb_ando(A,Y,In,t))
                    end
                    Y == B;
                cvx_end
                val2 = cvx_optval;

                cvx_clear
                cvx_quiet(true);
                cvx_begin
                    if cplx
                        variable X(n,n) hermitian
                        variable Y(n,n) hermitian
                    else
                        variable X(n,n) symmetric
                        variable Y(n,n) symmetric
                    end
                    if t >= 0 && t <= 1
                        maximize (lieb_ando(X,Y,In,t))
                    else
                        minimize (lieb_ando(X,Y,In,t))
                    end
                    X == A;
                    Y == B;
                cvx_end
                val12 = cvx_optval;

                QtAB = lieb_ando(A,B,In,t);
                [p,q] = rat(t);
                fprintf('n=%d, t=%d/%d, cplx=%d: ',n,p,q,cplx);

                % Thresholds are generous here to account for inaccurate
                % solving of SDPs (solver-dependent). In most cases,
                % solving the SDP should return a value that's within
                % 1e-8 or 1e-9 of true value (SDP representation in
                % lieb_ando is exact, there is no approximation, the error
                % can only come from solving the SDP)
                assert(abs(QtAB-val1) <= 1e-4,'Test failed lieb_ando n=%d, cplx=%d',n,cplx);
                assert(abs(QtAB-val2) <= 1e-4,'Test failed lieb_ando n=%d, cplx=%d',n,cplx);
                assert(abs(QtAB-val12) <= 1e-4,'Test failed lieb_ando n=%d, cplx=%d',n,cplx);
                fprintf('OK\n');
            end
        end
    end
    
    fprintf('lieb_ando test OK\n');

end

% Generate a random positive semidefinite matrix
function X = randPSD(n,iscplx)
    A = randn(n);
    if iscplx
        A = A +1i*randn(n);
    end
    X = A*A';
end

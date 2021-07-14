function cvxquad_tests()

    % CVXQUAD_TESTS
    % Run some unit tests on cvxquad's functions
    % These are very simple tests and are designed to quickly detect any
    % trivial problem with the code

    rng(0);

    test_matrix_geo_mean_hypo_cone();
    test_lieb_ando();
    test_op_rel_entr_epi_cone();
    test_trace_logm();
    test_quantum_entr();
    test_quantum_rel_entr();
    
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
            for mk = [1 3]
                for apx=[-1 0 1]
                    if mk == 1 && apx == 0
                        % When (m,k) = (1,1) we only check apx=\pm 1 and
                        % verify that the result is a valid upper/bound on
                        % the correct result. We do not check the quality
                        % of the result. For this setting we skip apx = 0
                        % (Pade) as this approximation is neither an upper
                        % nor a lower bound
                        continue;
                    end
                    cvx_quiet(true);
                    cvx_begin
                        if cplx
                            variable T(n,n) hermitian
                        else
                            variable T(n,n) symmetric
                        end
                        minimize trace(T)
                        {A,B,T} == op_rel_entr_epi_cone(n,cplx,mk,mk,eye(n),apx);
                    cvx_end
                    DopAB = A^(1/2)*logm(A^(1/2)*inv(B)*A^(1/2))*A^(1/2);
                    err = (T - DopAB)/norm(DopAB); % matrix
                    fprintf('n=%d, cplx=%d, (m,k)=(%d,%d), apx=%d, eig(err) in [%.4f,%.4f]: ',n,cplx,mk,mk,apx,min(eig(err)),max(eig(err)));
                    assert( min(eig(apx*err)) >= -1e-6,'Test failed (bound) op_rel_entr_epi_cone n=%d, cplx=%d, min(eig(apx*err))=%.4e',n,cplx,min(eig(apx*err)));
                    if mk >= 3
                        % Tolerance 1e-2 was set by inspection
                        assert(norm(err) <= 1e-2,'Test failed op_rel_entr_epi_cone n=%d, cplx=%d, error=%.4e',n,cplx,norm(err));
                    end
                    fprintf('OK\n');
                end
            end
        end
    end
    
    fprintf('op_rel_entr_epi_cone test OK\n');
end

function test_trace_logm()
    fprintf('---- TESTING trace_logm ----\n');
    nvec = [3 5 10];
    for n=nvec
        for cplx=[0 1]
            A = randPSD(n,cplx); A = A/trace(A);
            C = randPSD(n,cplx);
            for mk = [1 3]
                for apx=[-1 0 1]
                    if mk == 1 && apx == 0
                        % When (m,k) = (1,1) we only check apx=\pm 1 and
                        % verify that the result is a valid upper/bound on
                        % the correct result. We do not check the quality
                        % of the result. For this setting we skip apx = 0
                        % (Pade) as this approximation is neither an upper
                        % nor a lower bound
                        continue;
                    end
                    cvx_quiet(true);
                    cvx_begin
                        if cplx
                            variable X(n,n) hermitian
                        else
                            variable X(n,n) symmetric
                        end
                        maximize (trace_logm(X,C,mk,mk,apx))
                        X == A;
                    cvx_end
                    trlog = trace_logm(A,C);
                    relerr = (cvx_optval - trlog)/abs(trlog);
                    fprintf('n=%d, cplx=%d, (m,k)=(%d,%d), apx=%d, err=%.4f: ',n,cplx,mk,mk,apx,relerr);
                    % Check it is a true bound depending on value of apx
                    assert(apx*relerr >= -1e-6,'Test failed (bound) trace_logm n=%d, cplx=%d, apx=%d, error=%.4e',n,cplx,apx,relerr);
                    if mk >= 3
                        % Tolerance 1e-2 was set by inspection
                        assert(abs(relerr) <= 1e-2,'Test failed trace_logm n=%d, cplx=%d, apx=%d, error=%.4e',n,cplx,apx,relerr);
                    end
                    fprintf('OK\n');
                end
            end
        end
    end
end

function test_quantum_entr()
    fprintf('---- TESTING quantum_entr ----\n');
    nvec = [3 5 10];
    for n=nvec
        for cplx=[0 1]
            A = randPSD(n,cplx); A = A/trace(A);
            for mk = [1 3]
                for apx=[-1 0 1]
                    if mk == 1 && apx == 0
                        % When (m,k) = (1,1) we only check apx=\pm 1 and
                        % verify that the result is a valid upper/bound on
                        % the correct result. We do not check the quality
                        % of the result. For this setting we skip apx = 0
                        % (Pade) as this approximation is neither an upper
                        % nor a lower bound
                        continue;
                    end
                    cvx_quiet(true);
                    cvx_begin
                        if cplx
                            variable X(n,n) hermitian
                        else
                            variable X(n,n) symmetric
                        end
                        maximize (quantum_entr(X,mk,mk,apx))
                        X == A;
                    cvx_end
                    HA = quantum_entr(A);
                    relerr = (cvx_optval - HA)/abs(HA);
                    fprintf('n=%d, cplx=%d, (m,k)=(%d,%d), apx=%d, err=%.4f: ',n,cplx,mk,mk,apx,relerr);
                    % Check it is a true bound depending on value of apx
                    assert(apx*relerr >= -1e-6,'Test failed (bound) quantum_entr n=%d, cplx=%d, apx=%d, error=%.4e',n,cplx,apx,relerr);
                    if mk >= 3
                        % Tolerance 1e-2 was set by inspection
                        assert(abs(relerr) <= 1e-2,'Test failed quantum_entr n=%d, cplx=%d, apx=%d, error=%.4e',n,cplx,apx,relerr);
                    end
                    fprintf('OK\n');
                end
            end
        end
    end
end

function test_quantum_rel_entr()
    fprintf('---- TESTING quantum_rel_entr ----\n');
    for n=[2 3]
        for cplx=[0 1]
            A = randPSD(n,cplx); A = A/trace(A);
            B = randPSD(n,cplx); B = B/trace(B);
            for mk=[1 3]
                for apx = [-1 0 +1]
                    if mk == 1 && apx == 0
                        % When (m,k) = (1,1) we only check apx=\pm 1 and
                        % verify that the result is a valid upper/bound on
                        % the correct result. We do not check the quality
                        % of the result. For this setting we skip apx = 0
                        % (Pade) as this approximation is neither an upper
                        % nor a lower bound
                        continue;
                    end
                    
                    cvx_clear;
                    cvx_quiet(true);
                    cvx_begin
                        if cplx
                            variable X(n,n) hermitian
                        else
                            variable X(n,n) symmetric
                        end
                        minimize (quantum_rel_entr(X,B,mk,mk,apx))
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
                        minimize (quantum_rel_entr(A,Y,mk,mk,apx))
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
                        minimize (quantum_rel_entr(X,Y,mk,mk,apx))
                        X == A;
                        Y == B;
                    cvx_end
                    val12 = cvx_optval;

                    DAB = quantum_rel_entr(A,B);

                    % relative errors
                    err1 = (val1-DAB)/abs(DAB);
                    err2 = (val2-DAB)/abs(DAB);
                    err12 = (val12-DAB)/abs(DAB);

                    fprintf('n=%d, cplx=%d, (m,k)=(%d,%d), apx=%d, errors=(%.4f,%.4f,%.4f) ',n,cplx,mk,mk,apx,err1,err2,err12);

                    % Tolerance of 1e-2 was set by inspection
                    if mk >= 3
                        assert(abs(err1) <= 1e-2,'Test failed quantum_rel_entr n=%d, cplx=%d, apx=%d, error=%.4e',n,cplx,apx,abs(err1));
                        assert(abs(err2) <= 1e-2,'Test failed quantum_rel_entr n=%d, cplx=%d, apx=%d, error=%.4e',n,cplx,apx,abs(err2));
                        assert(abs(err12) <= 1e-2,'Test failed quantum_rel_entr n=%d, cplx=%d, apx=%d, error=%.4e',n,cplx,apx,abs(err12));
                    end

                    if apx ~= 0
                        % sign(err) should be the same as apx:
                        assert(apx*err1 >= -1e-6, 'Test failed (bound) quantum_rel_entr n=%d, cplx=%d, apx=%d, error=%.4e',n,cplx,apx,err1);
                        assert(apx*err2 >= -1e-6, 'Test failed (bound) quantum_rel_entr n=%d, cplx=%d, apx=%d, error=%.4e',n,cplx,apx,err2);
                        assert(apx*err12 >= -1e-6, 'Test failed (bound) quantum_rel_entr n=%d, cplx=%d, apx=%d, error=%.4e',n,cplx,apx,err12);
                    end

                    fprintf('OK\n');
                end
            end
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
        A = A + 1i*randn(n);
    end
    X = A*A';
end


function X = minSquareCholesky(M, RHS)
%% Solves the linear system by minimizing squares according to
%% Cholesky method
    M = full(M);
    opts.LT = true;
    L = chol(M'*M)';
    y = linsolve(L, M' * RHS, opts);
    X = linsolve(L', y);
end

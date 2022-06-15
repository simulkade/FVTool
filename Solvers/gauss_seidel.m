function X = gauss_seidel(M, RHS, X, opts)
% SYNOPSIS:
% gauss_seidel solves the linear system M.X = RHS with the iterative 
% Gauss-Seidel method. An optional relaxation factor can be used
%
%
% PARAMETERS:
% M : matrix of the linear system to solve
% RHS: right hand side
% solverSettings: struct containing options for the solver.
% includes the convergence criteria, the use of a relaxation factor, max and
% min number of iterations
%
% RETURNS:
% Solution X of the system
%
%
% EXAMPLE:
%
% SEE ALSO:
%

% Written by Lionel Sergent
% See the license file

arguments
    M double
    RHS double
    X double
    opts LSolverOpts = LSolverOpts()
end

sym_def_pos = true;
try
    chol(M);
catch
    sym_def_pos = false;
end

diag_dominant = true;
if any((2*abs(diag(M))) <= sum(abs(M),2))
    diag_dominant = false;
end

if ~diag_dominant && ~sym_def_pos && opts.warnings
    warning("Gauss-Seidel sufficient criterion not met")
end

X_prev = X;
N = numel(X);

    function cv = check_convergence()
        cv = true;
        % Little edge case to start
        if X == zeros(N,1) 
            if abs(M*X-RHS) == 0
                return
            else
                cv = false;
                return
            end
        end

        % Look at residues
        normed_residues = abs(M*X-RHS) ./ (abs(diag(M) .* X));
        if sum(normed_residues) > N * opts.rel_sum_residues_tol
            cv = false;
        end
        if max(normed_residues) > opts.rel_max_residues_tol
            cv = false;
        end

        % Look at changes
        normed_changes = 2 * abs(X - X_prev) ./ (abs(X) + abs(X_prev));
        if sum(normed_changes) > N * opts.rel_sum_change_tol 
            cv = false;
        end
        if max(normed_changes) > opts.rel_max_change_tol 
            cv = false;
        end
    end

iterations = 0;
iterations_since_cv = 0;
relax = opts.relaxation;

while (iterations_since_cv <  opts.post_iterations || ...
        iterations < opts.min_iterations)
    
    % One iteration = sweep Gauss-Seidel updates + 1 CV verification
    % Gauss-Seidel update
    for k = 1:opts.sweeps
    X_prev = X;
        for i = 1:N
            X(i) = relax/M(i,i) * (RHS(i) - M(i,:) * X) + X(i);
        end
    end

    % CV verification
    if check_convergence()
        iterations_since_cv = iterations_since_cv + 1;
    else
        iterations_since_cv = 0;
    end

    iterations = iterations + 1;

    if iterations >= opts.max_iterations
        msg = sprintf("Gauss-Seidel reached max iterations (%d) \n",...
            iterations);
        if ~sym_def_pos
            msg = msg + sprintf("Matrix is not symmetric definite \n");
        end
        if ~diag_dominant
            msg = msg + sprintf("Matrix is not diagonally dominant \n");
        end
        error(msg);
    end
end

end


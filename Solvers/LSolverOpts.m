classdef LSolverOpts
    properties
        solver_name = 'mldivide';
        warnings = true;

        % Properties for iterative solvers 
        initial_guess = 0;
        relaxation = 1; 
        sweeps = 10;  % number of sweeps per iteration
        max_iterations = 100;
        min_iterations = 3;
        post_iterations = 2;
        rel_sum_change_tol = 1e-2;
        rel_max_change_tol = 1e-1;
        rel_sum_residues_tol = 1e-4;
        rel_max_residues_tol = 1e-3;
    end

end


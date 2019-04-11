function err = vscmg_mt_sim_lin(K, params)
    %% Simulate case either with linearized or full dynamics

    %[~, A, B] = vscmg_mt_linearization(eye(9), eye(6), params);
    %[tsol, xsol] = ode15s(@(t, x) (A + B*K)*(x - params.x0), ...
    %                       [0, params.tf], params.x0);
    [tsol, xsol] = ode15s(@(t, x) vscmg_mt_dynamics(x, ...
                K*(x - params.x0), t, params), [0, params.tf], params.x0);
    
    %% Calculate error
    
    eps = xsol(:, 1:3);
    err = sum(trapz(tsol, eps.^2));
end


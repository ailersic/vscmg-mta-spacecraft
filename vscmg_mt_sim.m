function err = vscmg_mt_sim(K, params)
    %% Simulate case
    
    options = odeset('RelTol',1e-5);
    [tsol, xsol] = ode15s(@(t, x) vscmg_mt_dynamics(x, ...
        vscmg_mt_control(x, t, K, params), t, params), [0, params.tf], ...
        params.x0, options);

    %% Calculate error
    
    eps = xsol(:, 1:3);
    err = sum(trapz(tsol, eps.^2));
end

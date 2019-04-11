function convtime = vscmg_mt_sim_detumble(K, params)
    %% Simulate case
    
    options = odeset('RelTol',1e-5);
    [tsol, xsol] = ode15s(@(t, x) vscmg_mt_dynamics(x, ...
        vscmg_mt_control(x, t, K, params), t, params), [0, params.tf], ...
        params.x0, options);
    usol = vscmg_mt_control(xsol', tsol, K, params)';

    %% Calculate convergence time
    
    tol = 1e-3;
    isol = [];
    
    for i=1:length(tsol)
        if all(xsol(i, 1:6) < tol*ones(1, 6)) && all(usol(i, :) < ...
                tol*ones(1, 6))
            isol = i;
            break;
        end
    end
    
    if isempty(isol)
        convtime = params.tf;
    else
        convtime = tsol(isol);
    end
end

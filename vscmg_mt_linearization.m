function [K, A, B] = vscmg_mt_linearization(Q, R, params)
    %% Linearize system

    dx = 1e-6;
    du = 1e-6;

    A = zeros(length(params.x0), length(params.x0));
    B = zeros(length(params.x0), length(params.u0));

    for i=1:length(params.x0)
        dx0 = [zeros(i-1, 1); dx; zeros(length(params.x0)-i, 1)];
        A(:, i) = (vscmg_mt_dynamics(params.x0+dx0, params.u0, 0, ...
                                     params) - ...
                   vscmg_mt_dynamics(params.x0-dx0, params.u0, 0, ...
                                     params))/(2*dx);
    end

    for i=1:length(params.u0)
        du0 = [zeros(i-1, 1); du; zeros(length(params.u0)-i, 1)];
        B(:, i) = (vscmg_mt_dynamics(params.x0, params.u0+du0, 0, ...
                                     params) - ...
                   vscmg_mt_dynamics(params.x0, params.u0-du0, 0, ...
                                     params))/(2*du);
    end

    %% Check controllability
    
    if rank(ctrb(A, B)) == length(params.x0)
        disp('Controllable!')
    else
        disp('Not controllable!')
    end
    
    %% LQR Matrices and Riccati

    N = zeros(length(params.x0), length(params.u0));
    
    if min(eig([Q, N; N', R])) <= 0
        disp('Not solvable!')
    end

    K = -lqr(A, B, Q, R, N);
end
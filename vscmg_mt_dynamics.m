function [xdot, cdot, m, tau_mta, tau_vscmg, tau_ggs, tau_magres] ...
    = vscmg_mt_dynamics(x, u, t, params)
    %% Read input

    eps = x(1:3, :);
    eta = real(sqrt(1 - sum(eps.^2, 1)));
    
    q = [eta; eps];
    w = x(4:6, :);
    c = x(7:9, :);
    
    cdot = u(1:3, :);
    m = u(4:6, :);
    
    nt = size(x, 2);
    
    %% Physical constants
    
    mu = 6.67e-11*5.972e24; % m^3/s^2
    B0 = -7.9e15; % Wb m
    
    %% Orbital position and magnetic field
    
    n = sqrt(mu/(params.R^3)); % 1/s
    p_i = params.R*[cos(n*t'); ...
                    sin(n*t')*cos(params.incl); ...
                    sin(n*t')*sin(params.incl)];
    B_i = B0/(params.R^5)*[3*p_i(1,:).*p_i(3,:); ...
                           3*p_i(2,:).*p_i(3,:); ...
                           2*p_i(3,:).^2 - p_i(1,:).^2 - p_i(2,:).^2];
    
    C_bi = zeros(3, 3, nt);
    p_b = zeros(3, nt);
    B_b = zeros(3, nt);
    
    for i=1:nt
        C_bi(:, :, i) = [q(1, i)^2 + q(2, i)^2 - q(3, i)^2 - q(4, i)^2, ...
                         2*q(2, i)*q(3, i) - 2*q(1, i)*q(4, i), ...
                         2*q(2, i)*q(4, i) + 2*q(1, i)*q(3, i); ...
                         2*q(2, i)*q(3, i) + 2*q(1, i)*q(4, i), ...
                         q(1, i)^2 - q(2, i)^2 + q(3, i)^2 - q(4, i)^2, ...
                         2*q(3, i)*q(4, i) - 2*q(1, i)*q(2, i); ...
                         2*q(2, i)*q(4, i) - 2*q(1, i)*q(3, i), ...
                         2*q(3, i)*q(4, i) + 2*q(1, i)*q(2, i), ...
                         q(1, i)^2 - q(2, i)^2 - q(3, i)^2 + q(4, i)^2];
        
        p_b(:, i) = C_bi(:, :, i)*p_i(:, i);
        B_b(:, i) = C_bi(:, :, i)*B_i(:, i);
    end
    
    %% VSCMG
    
    G = zeros(3, 3, nt);
    L_vscmg = zeros(3, nt);
    
    for i=1:nt
        G(:, :, i) = params.Icmg*[cos(c(1, i))*cos(c(2, i))*c(3, i), ...
                                  -sin(c(1, i))*sin(c(2, i))*c(3, i), ...
                                  sin(c(1, i))*cos(c(2, i));
                                  0, ...
                                  cos(c(2, i))*c(3, i), ...
                                  sin(c(2, i)); ...
                                  -sin(c(1, i))*cos(c(2, i))*c(3, i), ...
                                  -cos(c(1, i))*sin(c(2, i))*c(3, i), ...
                                  cos(c(1, i))*cos(c(2, i))];
        
        L_vscmg(:, i) = params.Icmg*[sin(c(1, i))*cos(c(2, i))*c(3, i); ...
                                     sin(c(2, i))*c(3, i); ...
                                     cos(c(1, i))*cos(c(2, i))*c(3, i)];
    end
    
    %% Rigid body dynamics
    
    tau_mta = zeros(3, nt);
    tau_vscmg = zeros(3, nt);
    tau_ggs = zeros(3, nt);
    tau_magres = zeros(3, nt);
    tau = zeros(3, nt);
    wdot = zeros(3, nt);
    
    for i=1:nt
        tau_mta(:, i) = cross(m(:, i), B_b(:, i));
        tau_vscmg(:, i) = G(:, :, i)*cdot(:, i);
        tau_ggs(:, i) = 3*mu/(params.R^5)*cross(p_b(:, i), ...
                                                params.Isat*p_b(:, i));
        tau_magres(:, i) = cross(params.mres, B_b(:, i));

        tau(:, i) = tau_mta(:, i) + tau_vscmg(:, i) + ...
                    tau_ggs(:, i) + tau_magres(:, i);
        wdot(:, i) = params.invIsat*(tau(:, i) - ...
                                     cross(w(:, i), ...
                                     params.Isat*w(:, i) + ...
                                     L_vscmg(:, i)));
    end
    
    %% Assemble xdot
    
    qdot = zeros(4, nt);
    
    for i=1:nt
        qdotq = 0.5*quaternion([0; w(:, i)]')*quaternion(q(:, i)');
        [qdot(1, i), qdot(2, i), qdot(3, i), qdot(4, i)] = parts(qdotq);
    end
    
    xdot = [qdot(2:4, :); wdot; cdot];
end

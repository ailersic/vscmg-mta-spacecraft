function u = vscmg_mt_control(x, t, K, params)
    %% Read input

    eps = x(1:3, :);
    eta = real(sqrt(1 - sum(eps.^2, 1)));
    
    q = [eta; eps];
    w = x(4:6, :);
    c = x(7:9, :);
    
    nt = size(x, 2);
    
    %% Parse control inputs
    
    KcD = K(1);
    KmD = K(2);
    
    KcP = K(3);
    KmP = K(4);
    
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
    
    for i=1:nt
        G(:, :, i) = params.Icmg*[cos(c(1, i))*cos(c(2, i))*c(3, i), ...
                                  -sin(c(1, i))*sin(c(2, i))*c(3, i), ...
                                  sin(c(1, i))*cos(c(2, i)); ...
                                  0, ...
                                  cos(c(2, i))*c(3, i), ...
                                  sin(c(2, i)); ...
                                  -sin(c(1, i))*cos(c(2, i))*c(3, i), ...
                                  -cos(c(1, i))*sin(c(2, i))*c(3, i), ...
                                  cos(c(1, i))*cos(c(2, i))];
    end
    
    %% Calculate control inputs
    
    KcD = diag([KcD, ...
                KcD, ...
                KcD*params.x0(9)*params.cdotmax(3)/params.cdotmax(1)]);
    KcP = diag([KcP, ...
                KcP, ...
                KcP*params.x0(9)*params.cdotmax(3)/params.cdotmax(1)]);
    
    cdot = zeros(3, nt);
    m = zeros(3, nt);
    
    for i=1:nt
        cdot(:, i) = -KcD*G(:, :, i)'*w(:, i) - KcP*G(:, :, i)'*eps(:, i);
        m(:, i) = -KmD*cross(B_b(:, i), w(:, i)) - ...
                   KmP*cross(B_b(:, i), eps(:, i));
        for j=1:3
            cdot(j, i) = max(min(cdot(j, i), params.cdotmax(j)), ...
                             params.cdotmin(j));
            m(j, i) = max(min(m(j, i), params.mmax(j)), params.mmin(j));
        end
    end
    
    u = [cdot; m];
end


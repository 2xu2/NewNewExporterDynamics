function xiprime = icebergsim(policy,xi, xi_C, xi_E, rho_CC, rho_EE)
% Simulation for iceberg cost next period
% The input variables are policy function and iceberg cost for a period
% Draw the iceberg cost tomorrow from the markov chain
u = rand(1, 'double');
if xi == 100
    xiprime = xi_E;
elseif xi == xi_E
    if u > rho_EE
        xiprime = xi_C;
    else
        xiprime = xi_E;
    end
else
    if u > rho_CC
        xiprime = xi_E;
    else
        xiprime = xi_C;
    end
end
xiprime = policy * xiprime + (1 - policy) * 100;
end
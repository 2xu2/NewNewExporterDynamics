function Pi = profit(X, Qvalue, Evalue, theta, Cstar, alphan, alphak, w, r)
% Profit function for a plant in new exporter dynamics of Ruhl, Willis (2013)
% X is the export decision. X = 0 and 1
% Qvalue and Evalue are realizations of simulated AR(1) processes
% theta is the elasticity of substitution
% Cstar is the size of world aggregate demand relative to domestic demand
% alphan and alphak are parameters of the Copp-Douglas production function
% w is wage, r is the rent price of capital 
% State space is X, Epsilon and Q

M = (1 + X * (Qvalue^theta) * Cstar)^(1/theta) * Evalue;
n = (r/(M * (w * alphak/(r * alphan))^(alphak*(theta - 1)/theta - 1) * ...
    alphak * (theta - 1) /theta))^(1/(alphan * (theta - 1)/theta + alphak * (theta - 1)/theta - 1));
% or
% n = (w * (w*alphak/(r * alphan))^(-alphak*(theta - 1)/theta) * ...
%    theta/(alphan * M * (theta - 1))) ^ (theta/(alphan * (theta - 1) + alphak * (theta - 1) - theta))
% or
%n = ((r* alphan/(w*alphak))^(alphak*(theta - 1)/ theta) * w * theta/(alphan *...
%    M * (theta - 1)))^(theta/((alphan + alphak) * (theta - 1) -theta));
k = w * alphak * n/(r * alphan);

Pi = M * n^(alphan * (theta - 1)/theta) * k^(alphak * (theta - 1)/theta) - w * n - r * k;
end




function [Pi] = profit(X, xi, Evalue, Qvalue, theta, Cstar, alphan, alphak, w, r)
% Profit function for a plant in new exporter dynamics of Ruhl, Willis (2013)
% X is the export decision. X = 0 and 1
% Qvalue and Evalue are realizations of simulated AR(1) processes (vectors)
% theta is the elasticity of substitution
% Cstar is the size of world aggregate demand relative to domestic demand
% alphan and alphak are parameters of the Copp-Douglas production function
% w is wage, r is the rent price of capital 
% State space is X, Epsilon and Q
% xi is the iceberg cost. Note when X = 0, value of xi doesn't matter

% Not fastest, but clear
M = (1 + X * xi^(1 - theta) * Cstar .* (Qvalue.^theta)).^(1/theta) .* Evalue;
%define two cobb douglas coefficients
a = alphan * (theta - 1)/theta;
b = alphak * (theta - 1)/theta;
Pi = (1-a-b) * (a/w)^(a/(1-a-b)) * (b/r)^(b/(1-a-b)) .* M.^(1/(1-a-b));
end


